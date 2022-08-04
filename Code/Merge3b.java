import java.io.*;
import java.util.Arrays;
import java.util.Properties;

//THIS VERSION CAN BE RUN WITHOUT PARAMETERS, DEFAULTING TO ASSEMBLING ALL INDIVDUALS (1 TO 100000), WITH 5 THREADS

//April 13, 2018: VERSION 3 ALLOWS THE USER TO CHOOSE A LENGTH THRESHOLD, BELOW WHICH MERGED READS ARE NOT WRITTEN TO THE FILE. THIS IS TO REMOVE REPETITIVE READS THAT MAY BE ACCIDENTLALY MERGED

//13 Jan 2015: I modified the code to all for individual numbers >9999, but I did not change the program name.


//THIS VERSION DOES NOT FILTER OUT READS WITH N'S BUT INSTEAD TRIMMS OFF END N's AFTER MERGING...

//THE ASSEMBLE VERSION ADDITIONALLY DOES THE FOLLOWING:
//avoids use of param file
//autodetects read length (assumes that read length is constant)
//deletes original fastq files when finished merging
//filters out reads with one or more N's

////////////////////////////////////
//THIS VERSION TAKES A SET OF PAIRED READ FILES AND MERGES THEM, PRODUCING THREE OUTPUT FASTQ FILES, *M.FASTQ *U1.FASTQ *U2.FASTQ

public class Merge3b{

  public static void main(String[] args){
    try{

    	//default arguments
    	int begSample=1;
    	int endSample=100000;
    	int nThreads = 5;
    	int minMergedLengthAllowed = 50;

    	int nargs=args.length;

    	if(nargs>0){ begSample = Integer.parseInt(args[0]); }
    	if(nargs>1){ endSample = Integer.parseInt(args[1]); }
    	if(nargs>2){ nThreads = Integer.parseInt(args[2]); }
    	if(nargs>3){ minMergedLengthAllowed = Integer.parseInt(args[3]); }

		int nShortMergedReadsSkipped=0;
		
		//INITIALIZE VARIABLES WITH DEFAULT VALUES
		int maxlen=-1;		//what is the maximum length of the reads .. will be adjustd when read files are loaded
		int phredShift=33;	//assume phred shift of 33
		//int truncValue=999;
		int maxReads=999999999;
		double minpThresh=0.0000000001;
		boolean checkHeaders=false;
		int headerTrimAmt=0;
		boolean writeMeanQuals=false;

	  for(int sample=begSample; sample<=endSample; sample++){
		if(!new File("../"+getSampleName(sample)).exists()){continue;}
		String infolder = "../"+getSampleName(sample);
		String outstem = "../"+getSampleName(sample)+"/"+getSampleName(sample);

		System.out.println("Merging reads in folder "+infolder);

		int nReadsWithN=0;
		//CREATE A LIST OF THE R1 FASTQ FILES AND IDENTIFY WHICH ONES HAVE A MATE
		File dir = new File(infolder);

		//create filter for dir that contains only those files with .fastq and *_R1_*
		FilenameFilter filter = new FilenameFilter() {
		    public boolean accept(File dir, String name) {
		        return (name.endsWith(".fastq") && name.indexOf("_R1_")>0);
		    }
		};
		//File array to hold fastq files in dir as File objects
		File prelimFileNames[] = dir.listFiles(filter);
		boolean foundMate[] = new boolean[prelimFileNames.length];
		int nFilePairs=0;
		for(int i=0; i<prelimFileNames.length; i++){
			if(new File(infolder+"/"+prelimFileNames[i].getName().replace("_R1_","_R2_")).exists()){
				nFilePairs++;
				foundMate[i]=true;
			}else{
				System.out.println("Warning: Could not find a mate for :"+prelimFileNames[i].getName());
			}
		}
		
		if(nFilePairs==0){
			System.out.println("No valid paired read files were found. Skipping individual.");
			continue;
		}

		String inFileNamesR1[] = new String[nFilePairs];
		String inFileNamesR2[] = new String[nFilePairs];
		int counter=0;
		for(int i=0; i<nFilePairs; i++){
			if(foundMate[i]){
				inFileNamesR1[counter]=prelimFileNames[i].getPath();
				inFileNamesR2[counter]=inFileNamesR1[i].replace("_R1_","_R2_");
				System.out.println(inFileNamesR1[counter]+"\t"+inFileNamesR2[counter]);
				counter++;
			}
		}
		
		//identify the maximum read length (assume first read in each file is the same length as the remainder in that file)
		maxlen=-1;
		for(int i=0; i<nFilePairs; i++){
			BufferedReader brTemp = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(inFileNamesR1[i]) ) ));
			String temp=brTemp.readLine();
			temp=brTemp.readLine();
			if(temp.length()>maxlen){maxlen=temp.length();}
			brTemp.close();
		}		
				
		//CREATE OUTFILES
		BufferedWriter bwM = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(outstem+"_M.fastq") ) )); //outfile containing only merged reads
		BufferedWriter bwU1 = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(outstem+"_U1.fastq") ) )); //outfile containing only unmerged reads (read1)
		BufferedWriter bwU2 = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(outstem+"_U2.fastq") ) )); //outfile containing only unmerged reads (read2)
		

		//////////////////////////////////////
		//OPEN READ FILES AND CREATE OUTFILES
		//////////////////////////////////////
		int currFileIndex=0;
 		BufferedReader br1 = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(inFileNamesR1[currFileIndex]) ) ));	//forward read file
 		BufferedReader br2 = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(inFileNamesR2[currFileIndex]) ) ));	//reverse read file


		///////////////////////////
		//SETUP THE MULTITHREADING
		///////////////////////////
		MergeNode listBeg[] = new MergeNode[nThreads-1];
		MergeNode listEnd[] = new MergeNode[nThreads-1];
		int listLens[] = new int[nThreads-1];
		int currList=0;
		int maxListLen=20;

		Thread threads[] = new Thread[nThreads-1];

		int readNumb=0;
		int goodReadNumb=0;

		String tempS;
	
		MergeNode node = new MergeNode();
		node.head1=br1.readLine();
		node.head2=br2.readLine();
		readNumb++;
		
		while(node.head1.indexOf(":Y:")>0){
			tempS=br1.readLine();		//burn seqs
			tempS=br2.readLine();		//burn seqs
			tempS=br1.readLine();		//burn +
			tempS=br2.readLine();		//burn +
			tempS=br1.readLine();		//burn quals
			tempS=br2.readLine();		//burn quals
			
			node.head1=br1.readLine();	//get the next header
			node.head2=br2.readLine();	
			readNumb++;

			if(node.head1==null){break;}
		}

		int nLengths[] = new int[2*maxlen];	//maximum length is 2*len-1, 2*len th spot is for unmerged count	
		double pbinom[][] = new double[maxlen+1][maxlen+1]; //probability of k matches in n overlapping bases, where len is the maximum overlap
		System.out.println("Precalculating Probabilities...");
		//precompute the binomial probabilities for all possible cases:
		for(int n=1; n<=maxlen; n++){
			for(int k=0; k<=n; k++){
				for(int x=k; x<=n; x++){
					pbinom[n][k]+=Math.exp(choose(n,x,true) + Math.log(Math.pow(0.25,x)) + Math.log(Math.pow(1-0.25,n-x)));	//note that working in logs avoids overflow for factorial(>171)
				}	
			}
		}
		
		int nWritten=0;

		while(node.head1!=null){
			goodReadNumb++;
			//System.out.println("Evaluating Read "+goodReadNumb);
			node.read1=br1.readLine();	//read sequence
			node.read2=getRevComp(br2.readLine());
			tempS=br1.readLine();			//read in + line and ignore
			tempS=br2.readLine();
			node.qual1=br1.readLine();
			node.qual2=getRev(br2.readLine());

			if(readNumb%1000==0){
				System.out.print("\rEvaluating Read "+readNumb);
				for(int i=0; i<nThreads-1; i++){
					System.out.print("     "+listLens[i]);
				}
				System.gc();
			}
			if(readNumb>maxReads){break;}

			//check to make sure that reads line up
			if(checkHeaders){
				if(!node.head1.substring(0,node.head1.length()-headerTrimAmt).equals(node.head2.substring(0,node.head2.length()-headerTrimAmt))){
					System.out.println("ERROR: Read pair coordinates do not correspond:"+node.head1+"\t"+node.head2);
					return;
				}
			}

			//MASTER IS DONE READING IN THE DATA, NOW NEED TO PUT THE NODE IN THE LINKED LIST
			while(listLens[currList]>=maxListLen){
				//System.out.println("Main thread putting node in linked list");
				if(threads[currList]==null){
					//initialize the thread
					//System.out.println("Initializing thread "+currList);
					threads[currList] = new Thread(new MergeThread(listBeg[currList], pbinom, minpThresh,phredShift));
					//System.out.println("done.");
					threads[currList].start();
					//System.out.println("done starting run function.");
				}
				//write to file the processed node data
				while(listBeg[currList].processed){
					//System.out.println("Writing to file "+currList);
					//write to file
					
//System.out.println("A\t"+listBeg[currList].read1);					
//System.out.println("A\t"+listBeg[currList].qual1);					
					
					if(listBeg[currList].read1.indexOf("N")>=0 || listBeg[currList].read2.indexOf("N")>=0){
//System.out.println("B\t"+listBeg[currList].read1);					
//System.out.println("B\t"+listBeg[currList].qual1);
						//allows N's but trim if on end
						int firstGood=0;
						int lastGood=0;
						//left end of read 1
						for(firstGood=0; firstGood<listBeg[currList].read1.length(); firstGood++){
							if(listBeg[currList].read1.charAt(firstGood)!='N'){
								break;
							}
						}
						//right end of read 1
						for(lastGood=listBeg[currList].read1.length()-1; lastGood>=0; lastGood--){
							if(listBeg[currList].read1.charAt(lastGood)!='N'){
								break;
							}
						}
//System.out.println("C\t"+listBeg[currList].read1);											
//System.out.println("C\t"+listBeg[currList].qual1);
						if(lastGood-firstGood<=20){
							listBeg[currList].read1="";
							listBeg[currList].qual1="";
							listBeg[currList].read2="";
							listBeg[currList].qual2="";
							listBeg[currList]=listBeg[currList].next;
							listLens[currList]--;
							continue;
						}
						listBeg[currList].read1=listBeg[currList].read1.substring(firstGood,lastGood+1);
						listBeg[currList].qual1=listBeg[currList].qual1.substring(firstGood,lastGood+1);
//System.out.println("D\t"+listBeg[currList].read1);											
//System.out.println("D\t"+listBeg[currList].qual1);
						
						if(!listBeg[currList].merged){
//	System.out.println(listBeg[currList].qual2+"\n");

							//left end of read 2
							for(firstGood=0; firstGood<listBeg[currList].read2.length(); firstGood++){
								if(listBeg[currList].read2.charAt(firstGood)!='N'){
									break;
								}
							}
							//right end of read 2
							for(lastGood=listBeg[currList].read2.length()-1; lastGood>=0; lastGood--){
								if(listBeg[currList].read2.charAt(lastGood)!='N'){
									break;
								}
							}
							if(lastGood-firstGood<=20){
								listBeg[currList].read1="";
								listBeg[currList].qual1="";
								listBeg[currList].read2="";
								listBeg[currList].qual2="";
								listBeg[currList]=listBeg[currList].next;
								listLens[currList]--;
								continue;
							}
							listBeg[currList].read2=listBeg[currList].read2.substring(firstGood,lastGood+1);
							listBeg[currList].qual2=listBeg[currList].qual2.substring(firstGood,lastGood+1);
//		System.out.println(listBeg[currList].qual2+"\n");

						}
						
						nReadsWithN++;
					}
					
					
					if(listBeg[currList].read1.length()>=minMergedLengthAllowed){		
						nWritten++;
						if(listBeg[currList].merged){
	//System.out.println("E\t"+listBeg[currList].read1);											
	//System.out.println("E\t"+listBeg[currList].qual1);
							nLengths[listBeg[currList].read1.length()]++;
							bwM.write(listBeg[currList].head1+"\n");
							bwM.write(listBeg[currList].read1+"\n");
							bwM.write("+\n");
							bwM.write(listBeg[currList].qual1+"\n");
						}else{
	//System.out.println("F\t"+listBeg[currList].read1);											
	//System.out.println("F\t"+listBeg[currList].qual1);
							nLengths[nLengths.length-1]++;
							bwU1.write(listBeg[currList].head1+"\n");
							bwU1.write(listBeg[currList].read1+"\n");
							bwU1.write("+\n");
							bwU1.write(listBeg[currList].qual1+"\n");	
							bwU2.write(listBeg[currList].head2+"\n");
							bwU2.write(listBeg[currList].read2+"\n");
							bwU2.write("+\n");
							bwU2.write(listBeg[currList].qual2+"\n");
						}
					}else{
						nShortMergedReadsSkipped++;
					}
					listBeg[currList]=listBeg[currList].next;
					listLens[currList]--;
				}
			
				currList++;
				currList=currList%(nThreads-1);
			}
			if(listLens[currList]==0){	//first time
				listBeg[currList]=node;
				listEnd[currList]=node;
			}else{
				listEnd[currList].next=node;
				listEnd[currList]=node;
			}
			listLens[currList]++;

			//Get ready for the next round
			node = new MergeNode();

			node.head1=br1.readLine();
			node.head2=br2.readLine();
			readNumb++;

			while(node.head1!=null && node.head1.indexOf(":Y:")>0){
				tempS=br1.readLine();		//burn seqs
				tempS=br2.readLine();		//burn seqs
				tempS=br1.readLine();		//burn +
				tempS=br2.readLine();		//burn +
				tempS=br1.readLine();		//burn quals
				tempS=br2.readLine();		//burn quals
			
				node.head1=br1.readLine();	//get the next header
				node.head2=br2.readLine();	
				readNumb++;

				if(node.head1==null){break;}
			}
		
			if(node.head1==null && currFileIndex+1<inFileNamesR1.length && new File(inFileNamesR1[currFileIndex+1]).exists()){
				br1.close();
				br2.close();
				new File(inFileNamesR1[currFileIndex]).delete();
				new File(inFileNamesR2[currFileIndex]).delete();
				//load the next file
				currFileIndex++;
				System.out.println("\n"+currFileIndex+"\t"+inFileNamesR1[currFileIndex]+"\t"+inFileNamesR2[currFileIndex]);
	 			br1 = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(inFileNamesR1[currFileIndex]) ) ));	//forward read file
	 			br2 = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(inFileNamesR2[currFileIndex]) ) ));	//reverse read file
				bwM.flush();
				bwU1.flush();
				bwU2.flush();
				node.head1=br1.readLine();
				node.head2=br2.readLine();
				readNumb++;
			}
		}
		
		br1.close();
		br2.close();		
		new File(inFileNamesR1[currFileIndex]).delete();
		new File(inFileNamesR2[currFileIndex]).delete();

		//if(nWritten==0){System.out.println("NO GOOD READS FOUND, NEED TO KILL PROGRAM AND START WITH NEXT INDIVIDUAL;"); return;}

		boolean skipping=false;
		
		for(int i=0; i<nThreads-1; i++){

			if(listEnd[i]!=null){
				listEnd[i].done=true;
			}else{
				System.out.println("WARNING: NO GOOD READS WERE FOUND!");
				skipping=true;
			}
		}
		
		if(skipping){continue;}

		//at the end finish up the writing for the remaining nodes
		int nRemainingLists=nThreads-1;
		while(nRemainingLists>0){
			//write to file the processed node data
			while(listLens[currList]>0){
				//write to file
				if(listBeg[currList].merged){
					nLengths[listBeg[currList].read1.length()]++;
					bwM.write(listBeg[currList].head1+"\n");
					bwM.write(listBeg[currList].read1+"\n");
					bwM.write("+\n");
					bwM.write(listBeg[currList].qual1+"\n");
				}else{
					bwU1.write(listBeg[currList].head1+"\n");
					bwU1.write(listBeg[currList].read1+"\n");
					bwU1.write("+\n");
					bwU1.write(listBeg[currList].qual1+"\n");	
									
					bwU2.write(listBeg[currList].head2+"\n");
					bwU2.write(listBeg[currList].read2+"\n");
					bwU2.write("+\n");
					bwU2.write(listBeg[currList].qual2+"\n");
				}
				listBeg[currList]=listBeg[currList].next;
				listLens[currList]--;
				if(listLens[currList]==0){nRemainingLists--; break;}
			}
		
			currList++;
			currList=currList%(nThreads-1);
		}		

		bwM.flush();
		bwM.close();

		bwU1.flush();
		bwU1.close();

		bwU2.flush();
		bwU2.close();

	
		System.out.println();
		System.out.println(nReadsWithN+" reads contained one or more Ns (trimmed if on end)");
		
		System.out.println(nShortMergedReadsSkipped+"merged reads were determined to be too short and thus removed (not written in the output fastq files))");
		
		//System.out.println("Number of Pairs Pre-filter");
		//System.out.println("Total Number of Filtered Pairs = "+(countYes+countNo));
		//System.out.println("Number of overlapping pairs = "+countYes);
		//System.out.println("Number of non-overlapping pairs = "+countNo);

		//BufferedWriter bwX = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(outputFilename+"_I"+ind+"_lengths.txt") ) )); //merged read file	
		BufferedWriter bwX = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(outstem+"_lengths.txt") ) )); //merged read file
		for(int i=0; i<maxlen*2; i++){
			bwX.write(nLengths[i]+"\n");
		}
		bwX.flush();
		bwX.close();
	  }

    }catch(IOException ioe){System.out.println(ioe);}
  }

  static String getSampleName(int i){
  	String sampleName;
	if(i<10){sampleName="I000"+i;
	}else if(i<100){sampleName="I00"+i;
	}else if(i<1000){sampleName="I0"+i;
	}else{sampleName="I"+i;}
	//}else if(i<10000){sampleName="I"+i;
	//}else{System.out.println("Sample number exceeds 9999. Good Job cranking through the samples, but code not equipt!"); return "BAD";}  	
	return sampleName;
  }

  static String getRevComp(String seq){
  	String revComp="";
  	int len=seq.length();
  	for(int i=0; i<len; i++){
  		revComp+=seq.charAt(len-(i+1));
  	}
  	if(seq.length()!=revComp.length()){System.out.println("ERROR 1 IN getRevCOMP!!!");}
  
  	revComp=revComp.replace('A','Z');
  	revComp=revComp.replace('T','A');
  	revComp=revComp.replace('Z','T'); 
  	
  	revComp=revComp.replace('C','Q');
  	revComp=revComp.replace('G','C');
  	revComp=revComp.replace('Q','G');
  	
  	//REPLACE ALL OTHER BASES WITH 'N'
  	for(int i=0; i<len; i++){
  		char c=revComp.charAt(i);
  		if(c!='-' && c!='A' && c!='T' && c!='C' && c!='G'){
  			if(len>i+1){	revComp=revComp.substring(0,i)+'N'+revComp.substring(i+1);
  			}else{		revComp=revComp.substring(0,i)+'N'; }
		}
  	}  	
  	if(revComp.length()!=len){System.out.println("ERROR 2 IN getRevCOMP!!!");}
  		
	return revComp;
  }  

  static String getRev(String seq){
  	String rev="";
  	int len=seq.length();
  	for(int i=0; i<len; i++){
  		rev+=seq.charAt(len-(i+1));
  	}
 		
	return rev;
  }  

  static double choose(int n, int k){
  	return factorial(n)/(factorial(k)*factorial(n-k));
  }

  static double factorial(int x){
  	if(x<=1){	return 1.0;
  	}else{		return x*factorial(x-1); }
  }

  static double choose(int n, int k, boolean log){
	if(log){
  		return factorial(n,true)-((factorial(k,true)+factorial(n-k,true)));  	
  	}else{
  		return factorial(n)/(factorial(k)*factorial(n-k));
  	}
  }

  static double factorial(int x, boolean log){
  	if(log){
  		if(x<=1){	return Math.log(1.0);
  		}else{		return Math.log(x)+factorial(x-1, true); }
  	}else{
  		if(x<=1){	return 1.0;
  		}else{		return x*factorial(x-1, false); }  	
  	}
  }

}


class MergeNode{
	String head1,head2;
	String read1,read2;
	String qual1,qual2;
	boolean processed=false;
	boolean merged=false;
	boolean done=false;
	MergeNode next=null;
}

//Thread used as slave to merge reads
class MergeThread implements Runnable{
	Thread t;
	MergeNode node;
	int mini=-999;
	double minp=1;
	int miniB=-999;
	double minpB=1;
	int bestMatches=0;
	int matches=0;
	int len;
	double pmatches=0;
	double probF=-99;	//probability that the basecall is incorrect
	double probR=-99;
	double prF=-99;		//1-probF, or the probability that the base call is correct
	double prR=-99;
	int newPhred=-1;
	char newChar='!';
	double newProb=-999;
	int xF=-1;
	int xR=-1;
	double pbinom[][];
	double minpThresh;
	int phredShift;
	String newRead;
	String newQual;	
	
	public MergeThread(){}	//default constructor
	public MergeThread(MergeNode n, double pbinomX[][], double minpThreshX, int phredShiftX){
		pbinom=pbinomX;
		minpThresh=minpThreshX;
		phredShift=phredShiftX;
		t = new Thread(this);
		node=n;
	}
	public void run(){
		try{
			while(node.done==false){	 //node will self terminate when the master identifies that it is done reading/writing
				//System.out.println("Merging node "+node.head1);
				mini=-999;
				minp=1;
				miniB=-999;
				minpB=1;
				bestMatches=0;

				len=node.read1.length();
				/////////////////////////////////////////////////////
				//TRY OVERLAP
				/////////////////////////////////////////////////////
				for(int i=(len-1); i>=-(len-1); i--){
					matches=0;
					for(int j=Math.max(i,0); j<len+Math.min(0,i); j++){	
						if(node.read1.charAt(j)==node.read2.charAt(j-i)){matches++;}
					}
					pmatches=pbinom[Math.min(len-i,len+i)][matches];	

					if(pmatches<minp && i+len>0){
						minpB=minp;
						miniB=mini;
						minp=pmatches;
						mini=i;
						bestMatches=matches;
					}
					//COULD BE USED TO SPEED UP if(minp<minpThresh){break;}
				}

				int predFragLen=mini+len;

		        	//predicted to overlap
				if( minp<minpThresh){
					newRead="";
					newQual="";
					for(int x=0; x<len+mini; x++){

						xF=x;
						xR=x-mini;		

						if(xF>=0 && xF<=(len-1) && xR>=0 && xR<=(len-1)){ //consensus is being taken between the two corresponding base calls
							
							probF=(Math.pow(10.0,(-((int)(node.qual1.charAt(xF))-phredShift)/10.0)));	//p(!A) if A is base called
							probR=(Math.pow(10.0,(-((int)(node.qual2.charAt(xR))-phredShift)/10.0)));
							prF=1-probF;
							prR=1-probR;

							if(node.read1.charAt(xF)==node.read2.charAt(xR)){	//SAME BASE IS CALLED
								newChar=node.read1.charAt(xF);
								newProb= 1 - ( prF*prR )/( prF*prR + (0.25/0.75)*( 1.0 - prF - prR + prF*prR ) );
							}else{

								if(probR<=probF){
									newChar=node.read2.charAt(xR);
									newProb= 1 - ( prR*probF ) / ( prR*probF + (0.25/0.75)*(probR*(prR*probF+2.0)));
								}else if(probF<probR){
									newChar=node.read1.charAt(xF);
									newProb= 1 - ( prF*probR ) / ( prF*probR + (0.25/0.75)*(probF*(prF*probR+2.0)));
								}
							}					

						}else if(xF<0 || xF>(len-1)){ // only reverse strand is being used
							probR=(Math.pow(10.0,(-((int)(node.qual2.charAt(xR))-phredShift)/10.0)));
							newChar=node.read2.charAt(xR);
							newProb=probR;

						}else if(xR<0 || xR>(len-1)){ // only forward strand is being used
							probF=(Math.pow(10.0,(-((int)(node.qual1.charAt(xF))-phredShift)/10.0)));
							newChar=node.read1.charAt(xF);
							newProb=probF;

						}
						newRead+=newChar;
						newPhred=(int)(-Math.log10(newProb)*10);
						newQual=newQual+(char)(Math.min((int)(newPhred+phredShift),126));			//NOTE phredShift=33 CORRESPONDS TO SANGER CODE, note that a hard limit of asci 126 is imposed ('~')
					}

					node.merged=true;
					node.head1+="_M";
					node.read1=newRead;
					node.qual1=newQual;
					node.head2="";
					node.read2="";
					node.qual2="";
				}
				
				while(node.next==null){
					//System.out.println("run() waiting for more samples");
					t.sleep(10);
					if(node.done){node.processed=true; return;}
				}
				MergeNode tempNode=node;
				node=node.next;
				tempNode.processed=true;
			}
		}catch(InterruptedException e){System.out.println("Thread interrupted!");}
	}
}
