import scala.io.Source
import breeze.linalg._
import spark._

object Run_DFC_Spark{

  // Create breeze "CSCMatrix"
  private def createCSCFromFile(m: Int, n: Int, filename: String) = {
    val builder = new CSCMatrix.Builder[Double](m, n, 100339667)
    Source.fromFile(filename).getLines.map(_.split(' ')).map { elements =>
      (elements(0).toInt-1,elements(1).toInt-1,elements(2).toDouble)
    }.foreach { case(i,j,v) =>
      builder.add(i,j,v)
    }
    builder.result()
  }

  // Take a column slice of a CSCMatrix
  private def colSlice( A: CSCMatrix[Double], colIndex: Int ) : SparseVector[Double] = {
      val size = A.rows
      val rowStartIndex = A.colPtrs(colIndex)
      val rowEndIndex = A.colPtrs(colIndex + 1) - 1
      val capacity = rowEndIndex - rowStartIndex + 1
      val result = SparseVector.zeros[Double](size)
      result.reserve(capacity)
      var i = 0
      while( i < capacity ) {
         val thisindex = rowStartIndex + i
         val row = A.rowIndices(thisindex)
         val value = A.data(thisindex)
         result(row) = value
         i += 1
      }
      result
   }

  def main(args: Array[String]) {

    if (args.length < 9) {
      System.err.println("Usage: Run_DFC_Spark <master> <sparkHome> <jarfile> <trainfile> <testfile> <m> <n> <numsplits> <dataDir>")
      System.exit(1)
    }

    println("Master:           " + args(0))
    println("Spark home:       " + args(1))
    println("Jar file:         " + args(2))
    println("Train file:       " + args(3))
    println("Test file:        " + args(4))
    println("m:       " + args(5))
    println("n:        " + args(6))
    println("Number of splits: " + args(7))
    println("Data Dir:           " + args(8))

    val master = args(0)
    val sparkHome = args(1)
    val jarfile = args(2)
    val trainfile = args(3)
    val testfile = args(4)
    val m=args(5).toInt
    val n=args(6).toInt
    val numSplits = args(7).toInt
    val datadir = args(8)

    // time implementation
    val t_start = System.currentTimeMillis

    // fill training matrix
    val Mmatrix = createCSCFromFile(m, n, datadir+trainfile)
    println("Filled Mmatrix")

    // fill testing matrix
    val Lmatrix = createCSCFromFile(m, n, datadir+testfile)
    println("Filled Lmatrix")

    /* Partition matrix & run MF */
    val nparts = numSplits // set number of partitions
    val partsize = n/nparts // assume this is OK
    val sc = new SparkContext(master, "DFC_Spark")

    val M_array = (1 to Mmatrix.cols).map{ i => colSlice(Mmatrix,i-1)}
    val columnRdd = sc.parallelize(M_array, nparts)
    val acc = columnRdd.mapPartitionsWithSplit{ case(i, columns) =>
      val colArray = columns.toArray
      val builder = new CSCMatrix.Builder[Double](m, colArray.size)
      colArray.zipWithIndex.foreach { case(column, j) =>
        val myitr = column.array.iterator
        while(myitr.hasNext){
          val itr_val = myitr.next()
          builder.add(itr_val._1,j,itr_val._2)
        }
      }
      val localM = builder.result()
      Seq((i,sgd_base(localM))).toIterator
    }.collect().sortBy(_._1).map(_._2).toArray
    /* END MF partition */
    

    /* Recombine partitions with MF projection */
    println("Running MF projection")
    val L_proj = mf_proj(acc,n)
    val row_feats = L_proj._1
    val col_feats = L_proj._2
    /* END MF proj */


    /* Calculate error & return final results */
    var tr_err=0.0
    val tr_itr = Mmatrix.activeIterator
    while(tr_itr.hasNext){
      val curr_val = tr_itr.next()
      val predval_tr = (row_feats(::,curr_val._1._1):*col_feats(::,curr_val._1._2)).sum
      val curr_err = (predval_tr - curr_val._2)
      tr_err = tr_err + curr_err*curr_err
    }
    tr_err=math.sqrt(tr_err)

    var te_err=0.0
    val te_itr = Lmatrix.activeIterator
    while(te_itr.hasNext){
      val curr_val = te_itr.next()
      val predval_te = (row_feats(::,curr_val._1._1):*col_feats(::,curr_val._1._2)).sum
      val curr_err = (predval_te - curr_val._2)
      te_err = te_err + curr_err*curr_err
    }
    te_err=math.sqrt(te_err)
    /* END calc error */

    val t_end=System.currentTimeMillis

    // print final error
    println("DFC complete!")
    println("Training Error: " + tr_err/math.sqrt(Mmatrix.activeSize))
    println("Testing Error: " + te_err/math.sqrt(Lmatrix.activeSize))
    println("Elapsed Time: " + (t_end-t_start)/1000 + " seconds")
    /**********************************************************/
    /* END final results for DFC */

    sc.stop()
  }

  ////////// SGD BASE algorithm ////////////////
  def sgd_base(Mmatrix:CSCMatrix[Double]): (DenseMatrix[Double],DenseMatrix[Double]) = {
    
    /* Set parameters */
    val m = Mmatrix.rows
    val n = Mmatrix.cols
    val lrate = .04 // learning rate
    val k = .04 // parameter used to minimize over-fitting
    val min_impr = .0001 // min improvement required to continue iterating
    val init = .1 // initial value for features
    val rank = 10 // rank of feature vector
    val min_itrs=10

    /* Initialize */
    val minval = Mmatrix.min
    val maxval = Mmatrix.max
    val row_feats = DenseMatrix.zeros[Double](rank,m) += init
    val col_feats = DenseMatrix.zeros[Double](rank,n) += init
    var rmse = 2.0 // set rmse
    var rmse_prev = 2.0 // set previous rmse

    /* Find nonzero entries */
    val rowscols = Mmatrix.pairs.filter { case(k,v) => v != 0 }.iterator.map{case(k,v)=>k}.toArray
    val rows = rowscols.map{case(r,c)=>r}
    val cols = rowscols.map{case(r,c)=>c}
    val nonzeros = Mmatrix.pairs.filter { case(k,v) => v != 0 }.iterator.map{case(k,v)=>v}.toArray
    val num_nonzeros=nonzeros.size

    for(i<-0 until rank){
      println("rank: " + i)
      var t=0
      var impr=0.0
      while(t<min_itrs||impr>min_impr){
        var sq_err = 0.0
        for(j<-0 until num_nonzeros){
          // find predicted val
          val predvalvec = row_feats(::,rows(j)).t*col_feats(::,cols(j))
          var predval=predvalvec(0)
          // apply cut off
          if(predval<minval){
            predval=minval
          }
          if(predval>maxval){
            predval=maxval
          }
          // find error
          val err = nonzeros(j)-predval
          sq_err = sq_err + err*err
          // update row and col feats
          val new_row_feat = row_feats(i,rows(j))+lrate*(err*col_feats(i,cols(j))-k*row_feats(i,rows(j)))
          col_feats(i,cols(j))=col_feats(i,cols(j))+lrate*(err*row_feats(i,rows(j))-k*col_feats(i,cols(j)))
          row_feats(i,rows(j))=new_row_feat
        }
        // calculate rmse
        rmse_prev = rmse
        rmse = math.sqrt(sq_err/num_nonzeros)
        impr=rmse_prev-rmse
        t=t+1
      }
    }
    return (row_feats,col_feats)
  }

  ////////////////// MF projection ///////////////////
  def mf_proj(acc:Array[(DenseMatrix[Double],DenseMatrix[Double])], numcols:Int) : (DenseMatrix[Double],DenseMatrix[Double]) = {
    // project onto first partition ...
    val totalparts = acc.length // get the number of partitions
    val partsize = numcols/totalparts // assume okay
    val C = acc.head // get first element
    val X_A = C._1 // set output A to C.A
    val C_pinv_A = pinv((C._2).t)
    val C_pinv_B = pinv((C._1).t)
    val rank = X_A.rows
    val X_B = DenseMatrix.zeros[Double](rank,numcols)
    X_B(::,0 to partsize-1):=C._2 // set initial portion to first partition's B

    // iterate through the rest of the list ..
    acc.tail.zipWithIndex.foreach{ case(w,i) => 
      val M_hat_B = ((C._2*(C_pinv_A.t))*(C_pinv_B*(w._1).t))*w._2
      val curr_ind = (i+1)*partsize
      X_B(::, curr_ind to curr_ind+partsize-1):=M_hat_B
    }

    return (X_A,X_B)
  }
}
