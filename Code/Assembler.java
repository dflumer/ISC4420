import java.io.*;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Arrays;

//27 Jan 2015: Note, added feature to allow maxiumum number of reads to be considered. TO use this feature, create a file called maxReadsToConsider.txt with the number you want in a single line

//ARL 16 Oct 2013
//THIS ASSEMBLER BEGINS WITH A DIVERGENT MATCHING APPROACH USING SPACED KMERS FROM A SET OF REFERENCES
//KMERS FROM READS WITH DIVERGENT MATCHES ARE ADDED TO A LOCAL REFERENCE HASH TABLE USED FOR QUICK LOOKUP
//EACH KMER (OF EITHER TYPE) IS ALLOWED A MAXIMUM NUMBER OF MATCHES, THEN IS DISABLED
//ALSO ONCE A SPACED KMER GETS A VERIFIED MATCH, IT IS DISABLED
//THE SIZE OF THE SPACED KMER LIST SHOULD SHRINK WHILE THE SIZE OF THE LOCAL REFERENCE HASH TABLE GROWS
//REPETITIVE ELEMENTS ARE DEALT WITH BY IDENTIFYING THEM IN THE READS FILE AND NOT ALLOWING LOCAL KMERS CONTAINING THE ELEMENTS TO BE HASHED

//INPUT IS A FASTQ FILE, AND A SPACED KMER FILE (FROM IdentifySpacedKmers7.java)


//!*!*!*!*!* IDEAS TO TRY: 
//(1) check ends of read first, then middle...also try not checking each position (only odd or onlyeven...etc.) 
// measure success if keep trying after ends fail


public class Assembler{
	static Kmer spacedKmers=null;					//pointer to linked list containing spacedKmers
	static Kmer refKmers=null;						//pointer to linked list containing reference kmers
//	static Kmer localKmers=null;					//pointer to linked list containing local kmers

	static HashMap spacedMap = new HashMap();		//hash map containing the spaced kmers found in the *.kmers files...used on ensure each spaced kmers is ...
	static HashMap refMap = new HashMap();			//hash map containing exact match kmers from reference sequences
	static HashMap localMap = new HashMap();		//hash map containing exact match kmers from the species being assembled
	static HashMap contigPairMap = new HashMap();	//hash map containing contigID pairings, will be used to combine contigs
	static HashMap nmerMap = new HashMap();			//hash map containing all possible Nmers, will be used to shortcut eliminate spaced kmer lookups
	static HashMap contigCountMap = new HashMap();	//hash map with key=contigID, value= nReadsMappedInContig
	static HashMap relContigCountMap = new HashMap();	//hash map with key=contigID, value= nReadsMappedInContig/AvgReadsMapped

	static HashMap copyMap = new HashMap();	//used to count the number of copies of all 15mers observed in the reads...used to filter out repetitive elements on the fly

	static HashMap bad15merMap = new HashMap();		//just the copyMap members with more merCopies than the threshold

	static String refStems[];						//array of file stems for the different references
	static String refSeqs[][];						//array of reference sequences used for extensive match check
	static int refSeqLens[][];						//lengths of sequences in array
	static int siteToPos[][][];						//tables giving reference sequence positions for each alignment position... [loc][ref][site] ... -1 => NA
													//the tables were created so all corresponding reference positions could be obtained in order to determine the best matching reference during extensive (100bp) comparisons
	static int NmerKmerMatchesR[];					//for a given read, this tallies the number of times each unspacedSpaced kmer was found
	static int NmerKmerMatchesF[];					//for a given read, this tallies the number of times each unspacedSpaced kmer was found

	static Kmer kmerRefs[];							//array storing the kmer to which each read was matched...used to filter out reads matching before their kmer was disabled

	static int nUnspacedSpacedKmers=0;				//the number of 20 mers in kmer file that are not actually spaced
	static int nSpacedSpacedKmers=0;				//the number of 20 mers in kmer file that are not actually spaced
	static int nSpacedKmers=0;
	static int nRefKmers=0;
	static int nLocalKmers=0;
	
	static int nContigs=0;
	static int nContigPairs=0;
	
	static int localK=60;
	static int alignK=20;
	static int refK=20;
	
	static int N=8;									//small kmer used to avoid time-consuming spaced kmer search. E.g. if 0 of the 7mers within a given kmer exist in a read, then the kmer has only a 2% chance of matching to the read.
	
	static int kmerPosInRead=0;

	static int matchThreshold1= 17;					//number of matches out of 20 required for match to divergent reference using spaced kmers
	static double matchThreshold2= 0.55;			//Percent of matches required for final match check
	static int retireThreshold4=50;					//maximum depth for kmer chain

	static int kmerCoverageThreshold=3;				//in how many reads must a kmer be found to be used for local kmer mapping? ..!!this helps avoid assembles that derail due to rare chimeric sequences

	static int merCountThreshold=300;				//to be determined based on 15mer profile in reads file...used to prevent localKmer hashing when localKmer contains repetitive element
	static int contigAbundanceThreshold=10;			//contigs with 10x as many reads as the average are temporarily disabled (reads matching will not be added, but may be considered in the future if the other contigs catch up in abundance)

	static int readsAddedThreshold=1;

	static int mappedSumThreshold = 10;
	static double mappedRatioThreshold = 0.5;

	static int nMapped[][];
	static double nMappedRatio[][];
	static int nMappedSum[];
	
	static boolean skipLocus[];

	static int nLoci=0;
	static int nRefs=0;
	static int nSites=0;
	
	static int READID=0;
	
	public static void main(String[] args){
System.out.println("If you need to remove a locus, create the file LociToSkip.txt in the code folder and type in the locus number.");
	  try{
		long begTime = new Date().getTime();

		long tempTime;		
		
		String ind = args[0];	//e.g. 31252
		String readFileStem = "../I"+ind+"/I"+ind;	//e.g. ../I31252/I31252

		String refsFiles = "UNKNOWN";					//name of file containing the name stems of the reference files (3 files per reference, *.kmers, *.seqs, *.poss)

		//FIND REFERENCE LIST FILE
		File[] fileList = new File("../References").listFiles();
    	for(int i=0; i<fileList.length; i++){
    		String filename=(""+fileList[i]);
    		if(filename.endsWith(".txt")){	//first text file found in this folder is assumed to be the references file
    			refsFiles = filename;
    			break;
    		}
    	}
    	if(refsFiles.equals("UNKNOWN")){
    		System.out.println("ERROR: Could not find *.txt file in References folder.");
    		return;
    	}

    	//COUNT THE NUMBER OF LINES IN THE FILE
    	BufferedReader brZ = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(refsFiles) ) ));
    	String firstRefFile="";
    	nRefs=0;
    	String tempZ = brZ.readLine();
    	while(tempZ!=null && tempZ.length()>0){
    		if(firstRefFile.length()==0){
    			firstRefFile=tempZ;
    		}
    		nRefs++;
    		tempZ = brZ.readLine();
    	}
    	brZ.close();
    	System.out.println(nRefs+" references detected.");

    	//AUTODETECT THE NUMBER OF LOCI
		nLoci = 0;
		brZ = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(firstRefFile+".seqs") ) ));
    	tempZ = brZ.readLine();
    	while(tempZ!=null && tempZ.length()>0){
    		nLoci++;
    		tempZ = brZ.readLine();
    	}
    	brZ.close();
   		System.out.println(nLoci+" loci detected in references.");

		if(args.length>4){merCountThreshold=Integer.parseInt(args[4]);}	//optionally set the merCountThreshold to a different value after inspection of *.merCount file in R

		nMapped = new int[nRefs][nLoci];
		nMappedRatio = new double[nRefs][nLoci];
		nMappedSum = new int[nLoci];


		//declare reused variables
		int nLocalKmersNow=0;
		
		
		Kmer tempKmer;
		String tempSeq;
		String readF,readR,read;
		String qualF,qualR,qual;
		char orientation='F';

		int readLength=0;
		int kmerMatchesF=0;
		int kmerMatchesR=0;
		int kmerPosInReadF=-1;
		int kmerPosInReadR=-1;
		Read matchedReadF=null;
		Read matchedReadR=null;
		
		refStems = new String[nRefs];
		refSeqs = new String[nRefs][nLoci];
		refSeqLens = new int[nRefs][nLoci];

		//load in list of loci to skip, if available
		skipLocus = new boolean[nLoci];	//default false
		
		//this bit is used to skip over loci that are seen to cause assembly problems. if skipLocus is true then no spaced kmers will be hashed for the locus...
		BufferedReader brSKIP;
		if(new File("LociToSkip.txt").exists()){
			brSKIP = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("LociToSkip.txt") ) ));
			String tempSKIP = brSKIP.readLine();
			while(tempSKIP!=null){
				if(tempSKIP.length()!=0){
					skipLocus[Integer.parseInt(tempSKIP)-1]=true;
					System.out.println("\nNOTE: Skipping locus "+tempSKIP+" because LociToSkip.txt said to...");
				}
				tempSKIP=brSKIP.readLine();
			}
		}


		System.out.println("Loading information from references...");
		
		//SETUP THE SIXMERMAP
		setupNmerMap("");

		//LOAD THE REFERENCE SEQUENCES
		loadRefSeqs(refsFiles);
		
		//LOAD THE ALIGNMENT POSITIONS
		loadPoss(refStems);

		//LOAD KMER LIST
		loadKmers(refStems);

		NmerKmerMatchesF = new int[nUnspacedSpacedKmers];
		NmerKmerMatchesR = new int[nUnspacedSpacedKmers];

		//COUNT THE NUMBER OF LINES IN THE FILES
		String mer="";
		Integer merCount=null;
		String readFiles[] = {readFileStem+"_M.fastq",readFileStem+"_U1.fastq",readFileStem+"_U2.fastq"};
		int readsInFiles[] = new int[3];
		int totalReads=0;
		
		int maxMerCount=0;
		for(int R=0; R<3; R++){
			BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(readFiles[R]) ) ));
			String tempS=br.readLine();	//header

			while(tempS!=null){
				readsInFiles[R]++;
				totalReads++;

				 if(readsInFiles[R]%10000==0){
				 	System.out.print("\rProfiling reads for problematic 15mers to avoid..."+readsInFiles[R]);
				 }
				tempS=br.readLine(); //sequence
				
				if(totalReads<100000){						//only need to inspect first 100k reads to get an idea of which 15mers are super-abundant
					for(int i=0; i<tempS.length()-15; i++){
						mer=tempS.substring(i,i+15);
						merCount=(Integer)copyMap.get(mer);
						if(merCount==null){
							copyMap.put(mer,(Integer)1);
						}else{
							copyMap.put(mer,(Integer)(merCount+1));
							if(merCount+1>maxMerCount){
								maxMerCount=merCount+1;
							}
						}
					}
				}
				
				tempS=br.readLine(); //+
				tempS=br.readLine(); //quals
				tempS=br.readLine(); //next header
			}
			 System.out.println("\t"+readFiles[R]+" contains "+readsInFiles[R]+" reads.                                             ");
		}
		int totalReadsInFiles=readsInFiles[0]+readsInFiles[1]+readsInFiles[2];


		int printFreq1=Math.max(1,(int)Math.pow(10,Math.floor(Math.log10(totalReadsInFiles)-4)));
		int printFreq2=Math.max(1,(int)Math.pow(10,Math.floor(Math.log10(totalReadsInFiles)-3)));
		int updateFreq1=Math.max(1,(int)Math.pow(10,Math.floor(Math.log10(totalReadsInFiles)-2)));
		int updateFreq2=Math.max(1,(int)Math.pow(10,Math.floor(Math.log10(totalReadsInFiles)-1)));

		BufferedWriter bwMerCount = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(readFileStem+".merCount") ) ));	//store readID corresponding to the first occurance of each read sequence in file

		//merCountThreshold=100;
		//copy bad15Mers to new, leaner Hash Map
		Iterator it2 = copyMap.entrySet().iterator();
		while (it2.hasNext()) { 
			Map.Entry entry = (Map.Entry) it2.next(); 
			String key = (String)entry.getKey(); 
			int value = (Integer)entry.getValue();
			if(value>1){
				bwMerCount.write(key+"\t"+value+"\n");
			}
			if(value>merCountThreshold){
			//if(value>100){	//observed more than 100 in 100000 times
				bad15merMap.put(key,(Integer)value);
				bad15merMap.put(getRevComp(key),(Integer)value);	//put reverse complement in as well so don't have to check both orientations during hashLocalKmers
			}
		}
		bwMerCount.flush();
		bwMerCount.close();		

		// System.out.println("\t"+bad15merMap.size()+" 15mers will be disqualified based on count threshold of="+merCountThreshold+" per 100k reads.");		

		copyMap=null;	//remove large 15mer map to free up memory
		System.gc();


		//IDENTIFY REDUNDANT SEQUENCES TO AVOID

		boolean doneConsideringRead[] = new boolean[totalReadsInFiles];  ///true if mapped, or redundant				

		int dupCounts[] = new int[totalReadsInFiles];  ///0 if not yet seen (or if duplicate of a previously seen read),  >0 indicates number of times observed for first observed of its sequence			


	if(true){

		HashMap dupMap1 = new HashMap();			//uses first 20bp to identify some of the reads that are unique (if unique in first 20, than unique throughout) key=sequence.substring(0,20), value=readID
		HashMap dupMap2 = new HashMap();			//new version of dupmap uses read sequence itself to keep track of pcr duplicates (redundant sequences) key=sequence, value=readID

		Integer dupResult;
		String uniqueRC;
		int nDupsFound1=0;
		int nDupsFound2=0;
		int readsInspected=0;

		for(int R=0; R<3; R++){
			BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(readFiles[R]) ) ));
			String tempS,tempSRC;
			for(int i=0; i<readsInFiles[R]; i++){
				readsInspected++;
				 if(readsInspected%10000==0){
				 	System.out.print("\rIdentifying reads with redundant sequences (quick)..."+readsInspected);
				 }
				tempS=br.readLine();	//header
				tempS=br.readLine();	//read

				if(tempS.length()<=20){
					tempS=br.readLine();	//+
					tempS=br.readLine();	//qual
					continue;					
				}

				tempSRC=getRevComp(tempS);
				tempS=tempS.substring(0,20);
				tempSRC=tempSRC.substring(0,20);
			
				dupResult = (Integer)dupMap1.get(tempS);
				if(dupResult==null){
					//bwDup.write(readsInspected+"\n");
					dupMap1.put(tempS,readsInspected);
					dupMap1.put(tempSRC,readsInspected);
					dupCounts[readsInspected-1]=1;					
				}else{
					nDupsFound1++;
					dupCounts[dupResult-1]++;
				}
				
				tempS=br.readLine();	//+
				tempS=br.readLine();	//qual
			}
			br.close();
		}
		dupMap1=null;	//don't need dupMap anymore so make null to save RAM
		System.gc();
		 System.out.println("\n"+(int)(100*(readsInspected-nDupsFound1)/(double)readsInspected)+"% of reads were identified as unique using quick 20mer check... ");


		BufferedWriter bwDup = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(readFileStem+".dups") ) ));	//store readID corresponding to the first occurance of each read sequence in file

		readsInspected=0;
		for(int R=0; R<3; R++){
			BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(readFiles[R]) ) ));
			String tempS;
			for(int i=0; i<readsInFiles[R]; i++){
				readsInspected++;
				 if(readsInspected%10000==0){
				 	System.out.print("\rIdentifying reads with redundant sequences (exhaustive)..."+readsInspected);
				 }
				tempS=br.readLine();	//header
				tempS=br.readLine();	//read

				if(dupCounts[readsInspected-1]==1){
					tempS=br.readLine();	//+
					tempS=br.readLine();	//qual					
					doneConsideringRead[readsInspected-1]=false;
					continue;
				}

				dupResult = (Integer)dupMap2.get(tempS);
				if(dupResult==null){
					bwDup.write(readsInspected+"\n");
					dupMap2.put(tempS,readsInspected);
					dupMap2.put(getRevComp(tempS),readsInspected);
					dupCounts[readsInspected-1]=1;					
					doneConsideringRead[readsInspected-1]=false;
				}else{
					bwDup.write(dupResult+"\n");
					dupCounts[dupResult-1]++;
					nDupsFound2++;
					doneConsideringRead[readsInspected-1]=true;
				}

				tempS=br.readLine();	//+
				tempS=br.readLine();	//qual
			}
			br.close();
		}
		bwDup.flush();
		bwDup.close();
		dupMap2=null;	//don't need dupMap anymore so make null to save RAM
		System.gc();
		 System.out.println("\n"+(int)(100*(readsInspected-nDupsFound2)/(double)readsInspected)+"% of reads were identified as unique using quick and exhaustive checks... ");

		 System.out.println("\n"+(int)(100*nDupsFound2/readsInspected)+"% of the read sequences are redundant.");

	}


		//write down the number of reads in each file...this will be used later for allele phasing to connect polymorphisms spanned by unmerged read pairs...
		BufferedWriter bwRIF = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(readFileStem+"_readsInFiles.txt") ) ));
		bwRIF.write(readFiles[0]+"\t"+readsInFiles[0]+"\n");
		bwRIF.write(readFiles[1]+"\t"+readsInFiles[1]+"\n");
		bwRIF.write(readFiles[2]+"\t"+readsInFiles[2]+"\n");
		bwRIF.flush();
		bwRIF.close();
		
		 System.out.println("\t"+((new Date().getTime()-begTime)/1000)+" seconds elapsed so far.");

		Kmer kmerRefs[] = new Kmer[totalReadsInFiles];		

		int startReadID[] = new int[3];	//to keep track of which read the U1 and U2 begin on so contigs can be matched up later... also will be used to disable reads once they are matched...


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////// P R E L E M E N A R Y   A S S E M B L Y ////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		int nReadsAdded=999999;
		int totalReadsMapped=0;
		int round=0;
		
		int spacedKmerByPassedPercent=0;
		int localKmerByPassedPercent=0;
		
		BufferedWriter bw = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(readFileStem+".ass") ) ));	//preliminary assembly file

		//MAP ALL OF THE READS
		System.out.println("\nMapping reads...");
		while(nReadsAdded>readsAddedThreshold){
			round++;
			nReadsAdded=0;
			READID=0;
			int nReadsConsidered=0;	//accounts for duplicates
			for(int R=0; R<3; R++){
				startReadID[R]=READID;
				BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(readFiles[R]) ) ));	//ONLY UTILIZES ONE READ FILE
				String dummy=br.readLine();		//get the first header
				while(dummy!=null){
					READID++;
					nReadsConsidered+=dupCounts[READID-1];
					 if(round==1 && READID%printFreq1==0 || round>1 && READID%printFreq2==0){
					 	System.out.print("\r     Round "+round+", read "+READID+"  "+(int)(100*READID/(double)totalReadsInFiles)+"% of reads processed, ");
					 	if(round==1){
					 		System.out.print((int)(100*totalReadsMapped/(double)nReadsConsidered)+"% on-target ("+totalReadsMapped+" reads mapped), "+spacedKmerByPassedPercent+ "% spaced kmers bypassed.");
					 	}else{
					 		System.out.print((int)(100*totalReadsMapped/(double)totalReadsInFiles)+"% on-target ("+totalReadsMapped+" reads mapped), "+localKmerByPassedPercent+ "% local kmers bypassed.");							
					 	}
					 }
					
					if(READID%updateFreq1==0){
						if(round==1){
							//update skip list
							for(int i=0; i<nLoci; i++){
								nMappedSum[i]=0;
								double max=0;
								for(int j=0; j<nRefs; j++){
									nMappedSum[i]+=nMapped[j][i];
									if(nMapped[j][i]>max){max=nMapped[j][i];}
								}
								for(int j=0; j<nRefs; j++){
									if(max==0){
										nMappedRatio[j][i]=1;
									}else{
										nMappedRatio[j][i]=nMapped[j][i]/max;
									}
								}
							}
							//disabled kmers not useful
							tempKmer=spacedKmers;
							int nKmers=0;
							int nNotUseful=0;
							while(tempKmer!=null){
								nKmers++;
								//check to see if this kmer is still useful
								tempKmer.notUseful=nMappedSum[tempKmer.locus-1]>mappedSumThreshold && nMappedRatio[tempKmer.ref-1][tempKmer.locus-1]<mappedRatioThreshold;
								if(tempKmer.notUseful){nNotUseful++;}
								tempKmer=tempKmer.next;
							}
							spacedKmerByPassedPercent=(int)(100*nNotUseful/nKmers);

						}else if(READID%updateFreq2==0){
							Iterator it = localMap.entrySet().iterator();
							HashMap disableMap = new HashMap(); 
							while (it.hasNext()) { 
								Map.Entry entry = (Map.Entry) it.next(); 
								String key = (String)entry.getKey(); 
								Kmer value = (Kmer)entry.getValue();
								if(round-value.roundLastMatched>1 || (round-value.roundLastMatched==1 && value.readLastMatched<READID)){
									//need to disable this local kmer, because it can no longer be useful
									disableMap.put(key,null);
								}
							}

							it = disableMap.entrySet().iterator();
							while(it.hasNext()){
								Map.Entry entry = (Map.Entry) it.next();
								String key = (String)entry.getKey();
								localMap.remove(key);
							}
							nLocalKmersNow=localMap.size();
							localKmerByPassedPercent=(int)(100*(nLocalKmers-nLocalKmersNow)/nLocalKmers);
						}
						
						
						//write the coverage as a function of locus and reference
						BufferedWriter bwZ = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(readFileStem+".cov") ) ));	//unique readsMapped refxlocus
						for(int i=0; i<nLoci; i++){
							for(int j=0; j<nRefs; j++){
								bwZ.write(nMapped[j][i]+"\t");
							}
							bwZ.write("\n");
						}
						bwZ.flush();
						bwZ.close();

						//write coverage as a function of contig (used to detect issues with repetetive contig)
						BufferedWriter bwConCount = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(readFileStem+".con") ) ));	//number of reads mapped to each contig so far
						Iterator it = contigCountMap.entrySet().iterator();
						double avgContigCount=0;
					 	double nContigs=0;
						while (it.hasNext()) { 
							Map.Entry entry = (Map.Entry) it.next(); 
							Integer key = (Integer)entry.getKey(); 
							Integer value = (Integer)entry.getValue();
							if(value>2){	//ignore contigs with only one or two reads
								nContigs++;
								avgContigCount+=value;
							}
						}						
						avgContigCount/=nContigs;
						it = contigCountMap.entrySet().iterator();
						while (it.hasNext()) {
							Map.Entry entry = (Map.Entry) it.next(); 
							Integer key = (Integer)entry.getKey(); 
							Integer value = (Integer)entry.getValue();
							relContigCountMap.put(key,value/avgContigCount);	//number of reads in contig relative to the average (skipping 0 or 1)
							bwConCount.write(key+"\t"+value+"\t"+relContigCountMap.get(key)+"\n");
						}
						bwConCount.flush();
						bwConCount.close();
					}

					if(doneConsideringRead[READID-1]){
						dummy=br.readLine();	//read
						dummy=br.readLine();	//+
						dummy=br.readLine();	//quals
						dummy=br.readLine();	//next header
						continue;
					}

					readF=br.readLine();		//the read sequence

					readR=getRevComp(readF);	//the read sequence RevComp
					dummy=br.readLine();		//+

					qualF=br.readLine();		//quals
					qualR=getRev(qualF);		//quals

					readLength=readF.length();

					matchedReadF=null;	//these will be initialized when a potential match is made
					matchedReadR=null;

					//FIRST, CHECK THE LOCAL KMER HASH MAP (VERY FAST, NOT TOLERANT OF MISMATCHES)
					tempTime = new Date().getTime();
					matchedReadF=getKmer(readF,localMap,localK);
					matchedReadR=getKmer(readR,localMap,localK);

					if(matchedReadF!=null){
						matchedReadF.matchStatus=1;
						matchedReadF.kmer.roundLastMatched=round;
						matchedReadF.kmer.readLastMatched=READID;
					}
					if(matchedReadR!=null){
						matchedReadR.matchStatus=1;
						matchedReadR.kmer.roundLastMatched=round;	//note the match so this kmer will not be disabled soon
						matchedReadR.kmer.readLastMatched=READID;	//note the match so this kmer will not be disabled soon						
					}

					//SECOND,CHECK THE REFERENCE KMER HASH MAP (VERY FAST, NOT TOLERANT OF MISMATCHES)
					if(round==1 && matchedReadF==null && matchedReadR==null){
						tempTime = new Date().getTime();				
						matchedReadF=getKmer(readF,refMap,refK);
						if(matchedReadF!=null){					//now verify
							getMatches(readF, matchedReadF);
							if(matchedReadF.matchScore>matchThreshold2){matchedReadF.matchStatus=2;}else{matchedReadF=null;}
						}

						matchedReadR=getKmer(readR,refMap,refK);				
						if(matchedReadR!=null){												//now verify
							getMatches(readR, matchedReadR);
							if(matchedReadR.matchScore>matchThreshold2){matchedReadR.matchStatus=2;}else{matchedReadR=null;}
						}
					}

					//NEITHER IS VERIFIED YET, SO NEED TO CHECK THE SPACED KMER LIST (VERY SLOW, TOLERANT OF MISMATCHES)
					if(round==1 && matchedReadF==null && matchedReadR==null){			//skip exhaustive check if match was found in either orientation...

						//ATTEMPT TO RULE OUT KMER USING NMER MATCHES
						tempTime = new Date().getTime();
						Arrays.fill(NmerKmerMatchesF,0);
						Arrays.fill(NmerKmerMatchesR,0);
						int maxNmerKmerMatchesF=-1;
						int maxNmerKmerMatchesR=-1;

						for(int x=0; x<readLength-N; x++){
							Nmer nmer = (Nmer)nmerMap.get(readF.substring(x,x+N));
							if(nmer==null){continue;}
							nmer=nmer.next;		//skip the head node
							while(nmer!=null){
								NmerKmerMatchesF[nmer.kmerID]++;

								if(NmerKmerMatchesF[nmer.kmerID]>maxNmerKmerMatchesF){
									maxNmerKmerMatchesF=NmerKmerMatchesF[nmer.kmerID];
								}
								nmer=nmer.next;
							}

							nmer = (Nmer)nmerMap.get(readR.substring(x,x+N));
							if(nmer==null){continue;}
							nmer=nmer.next;		//skip the head node
							while(nmer!=null){
								NmerKmerMatchesR[nmer.kmerID]++;
								if(NmerKmerMatchesR[nmer.kmerID]>maxNmerKmerMatchesR){
									maxNmerKmerMatchesR=NmerKmerMatchesR[nmer.kmerID];
								}
								nmer=nmer.next;
							}
						}

						tempTime = new Date().getTime();
						//first check forward
						matchedReadF=getKmer(readF,NmerKmerMatchesF);
						if(matchedReadF!=null){
							getMatches(readF, matchedReadF);								

							if(matchedReadF.matchScore>matchThreshold2){matchedReadF.matchStatus=3;}else{ matchedReadF=null;}

						}
						//second check reverse
						matchedReadR=getKmer(readR,NmerKmerMatchesR);
						if(matchedReadR!=null){
							getMatches(readR, matchedReadR);

							if(matchedReadR.matchScore>matchThreshold2){matchedReadR.matchStatus=3;}else{ matchedReadR=null;}

						}
					}
			
					Read matchedRead=null;
					read=null;
					qual=null;

					//DETERMINE ORIENTATION OF BEST MATCH (IF ANY)
					if(matchedReadF==null && matchedReadR==null){
						//neither direction matched, so do nothing
					}else if(matchedReadF!=null && matchedReadR!=null){
						//both directions matched choose better
						if(matchedReadF.matchScore>=matchedReadR.matchScore){	//note that if matchStatus == 0 for both, then forward will be chosen by default
							matchedRead=matchedReadF;
							matchedRead.orientation='F';
							read=readF;
							qual=qualF;
						}else{
							matchedRead=matchedReadR;
							matchedRead.orientation='R';		
							read=readR;
							qual=qualR;	
						}
					}else if(matchedReadF!=null){
						matchedRead=matchedReadF;
						matchedRead.orientation='F';		
						read=readF;
						qual=qualF;	
					}else{
						matchedRead=matchedReadR;
						matchedRead.orientation='R';				
						read=readR;
						qual=qualR;
					}
			
			
					if(matchedRead!=null && matchedRead.contigID>=0){
						Integer result = (Integer)contigCountMap.get(matchedRead.contigID);
						if(result==null){
							contigCountMap.put(matchedRead.contigID,1);
						}else{
							contigCountMap.put(matchedRead.contigID,result+1);	
						}
			
						//disable the read if it maps to a contig with too high of abundance
						Double resultD = (Double)relContigCountMap.get(matchedRead.contigID);
						if(resultD!=null && resultD>contigAbundanceThreshold){
							matchedRead=null;					
						}
					}
			
					//DETERMINE WHETHER THE READ CONTAINS KMERS FROM MULTIPLE LOCI, IF SO kILL READ BEFORE IT WREAKS HAVOC
					String kmerSeq="";
					Kmer currKmer;
					if(matchedRead!=null && read!=null){
						for(int i=0; i<=read.length()-localK; i++){
							kmerSeq=read.substring(i,i+localK);
							currKmer = (Kmer)localMap.get(kmerSeq);
							if(currKmer!=null && currKmer.locus!=matchedRead.locusID){		//all reads with kmers belonging to more than one locus are disqualified here
								matchedRead=null;
								break;
							}
						}
					}
			
					if(matchedRead==null){
//						matchStatusCounter[0]++;	
					}else{						
						nReadsAdded+=dupCounts[matchedRead.readID-1];
						totalReadsMapped+=dupCounts[matchedRead.readID-1];
						
						//run along read sequence and hash the kmers  these are LOCAL KMERS
 
						if(matchedRead.kmer.depth<=retireThreshold4){
							hashLocalKmers(read,matchedRead,round);
						}

						kmerRefs[matchedRead.readID-1]=matchedRead.kmer;	//store kmer link so it can be used later to ignore reads descending form disabled kmers
						
						if(matchedRead.kmer.disabled || (matchedRead.kmer.parent!=null && matchedRead.kmer.parent.disabled)){
							System.out.println("\nRead from disabled kmer is going to be written!");
						}
						
						//note that the read is now mapped
						doneConsideringRead[READID-1]=true;
						
						nMapped[matchedRead.refID-1][matchedRead.locusID-1]++;

						bw.write(">loc"+matchedRead.locusID+"_read"+matchedRead.readID+"_ref"+matchedRead.refID+"_con"+matchedRead.contigID+"_ori"+matchedRead.orientation+"_aPos"+matchedRead.alignPos+"_"+"\n");
						bw.write(read+"\n");
						bw.write(qual+"\n");

						bw.flush();

					}
					dummy=br.readLine();		//get the next header

				}
				
			}
			bw.flush();
			// System.out.println("   "+nReadsAdded+" reads added in round "+round);
		}
		bw.flush();
		bw.close();

		//bwMER.flush();
		//bwMER.close();

		 System.out.println("\n\t"+((new Date().getTime()-begTime)/1000)+" seconds elapsed so far.");



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////// S U P E R C O N T I G   C R E A T I O N ////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////IDENTIFY CONTIG PAIRS USING READS IN *.ASS FILE///////////////////////
for(int loc=1; loc<=nLoci; loc++){
	System.out.print("\rIdentifying Supercontigs for locus "+loc+"...     ");
	
	//count the number of contigs for this locus...
	BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(readFileStem+".ass") ) ));
	String tempS=br.readLine();
	int nContigsThisLocus=0;
	HashMap contigHash= new HashMap();
	int contigID1=-1;
	while(tempS!=null){
		if(tempS.startsWith(">loc"+loc+"_")){	//e.g. >loc1_read10_ref2_con1_oriF_aPos-148_
			contigID1=Integer.parseInt(tempS.split("_")[3].substring(3));
			if(contigHash.get(contigID1)==null){
				nContigsThisLocus++;
				contigHash.put((Integer)contigID1,nContigsThisLocus-1);	//value corresponds to the array index below in contigPairs and contigIDs
			}
		}
		tempS=br.readLine();
	}
	br.close();
	 //System.out.print(""+nContigsThisLocus+" contigs found.\n");
	//now build contig pair matrix
	int contigPairs[][] = new int[nContigsThisLocus][nContigsThisLocus]; //default 0

	br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(readFileStem+".ass") ) ));
	tempS=br.readLine();

	contigID1=-1;
	String kmer="";
	Integer contigID2=-1;

	HashMap kmerHash = new HashMap();
	int contigIDs[] = new int[nContigsThisLocus];
	while(tempS!=null){
		if(tempS.startsWith(">loc"+loc+"_")){	//e.g. >loc1_read10_ref2_con1_oriF_aPos-148_
			contigID1=Integer.parseInt(tempS.split("_")[3].substring(3));
			contigIDs[(Integer)contigHash.get(contigID1)]=contigID1;	//need a way to retrieve the contigIDs from the matrix indexes below
			tempS=br.readLine();	//get sequence
			for(int i=0; i<=tempS.length()-localK; i++){	//look at every kmer in sequence
				kmer=tempS.substring(i,i+localK);
				if(kmer.indexOf("N")>=0){continue;} //skip kmers with N's
				contigID2=(Integer)kmerHash.get(kmer);

				if(contigID2==null){
					kmerHash.put(kmer,(Integer)contigID1);
				}else if(contigID1!=contigID2){		//found two contigs containing the same kmer, so connect
					contigPairs[(Integer)contigHash.get(contigID1)][(Integer)contigHash.get(contigID2)]++;
					contigPairs[(Integer)contigHash.get(contigID2)][(Integer)contigHash.get(contigID1)]++;
				} 
			}
		}
		tempS=br.readLine();
	}
	br.close();
	//if(loc>1){continue;}
	//print the current contig matrix
	for(int i=0; i<nContigsThisLocus; i++){
		for(int j=i+1; j<nContigsThisLocus; j++){
			if(contigPairs[i][j]>0){
				if(contigIDs[i]<contigIDs[j]){
					contigPairMap.put(contigIDs[i]+"_"+contigIDs[j],"-1_-1_-1_-1");
				}else{
					contigPairMap.put(contigIDs[j]+"_"+contigIDs[j],"-1_-1_-1_-1");				
				}
			}
		}
	}
}
System.out.println();		
		//COMBINE THE CONTIGS INTO SUPERCONTIGS	
		//create contigs and link them to the contig array
		Contig contigs[] = new Contig[nContigs];					//this should be equal to the number of reads that are matched directly to a non-local reference
		for(int i=0; i<nContigs; i++){
			contigs[i]=new Contig();
			contigs[i].contigID=i+1;
		}
		//create superContigs and link them to the superContig array
		SuperContig superContigs[] = new SuperContig[nContigs];	    //initially the number of contigs


		for(int i=0; i<nContigs; i++){ //nContigPairs; i++){
			superContigs[i]=new SuperContig();
			superContigs[i].superContigID=i+1;			
		}

		//assign sc ids
		Iterator it = contigPairMap.entrySet().iterator();
		int nSuperContigs=0;
		while (it.hasNext()) { 
			Map.Entry entry = (Map.Entry) it.next(); 
			String key = (String)entry.getKey(); 
			String value = (String)entry.getValue();
			Contig contigA=contigs[Integer.parseInt(key.split("_")[0])-1];
			Contig contigB=contigs[Integer.parseInt(key.split("_")[1])-1];

			if(contigA.superContig==null && contigB.superContig==null){	//create a new supercontig
				superContigs[nSuperContigs].firstContig=contigA;
				superContigs[nSuperContigs].lastContig=contigB;
				contigA.superContig=superContigs[nSuperContigs];
				contigB.superContig=superContigs[nSuperContigs];
				contigA.next=contigB;

				superContigs[nSuperContigs].locusID=Integer.parseInt(value.split("_")[3]);

				nSuperContigs++;				
			}else if(contigA.superContig==null){						//add contig A to the end of B's supercontig list
				contigA.superContig=contigB.superContig;
				contigA.superContig.lastContig.next=contigA;
				contigA.superContig.lastContig=contigA;
			}else if(contigB.superContig==null){						//add contig B to the end of A's supercontig list
				contigB.superContig=contigA.superContig;
				contigB.superContig.lastContig.next=contigB; 
				contigB.superContig.lastContig=contigB;				
			}else if (contigA.superContig!=contigB.superContig){	//both already have supercontigs, but they are not the same so need to collapse to supercontigs into one
				if(contigA.superContig.locusID!=contigB.superContig.locusID){
					System.out.println("WARNING: Two contings from different loci ("+contigA.superContig.locusID+","+contigB.superContig.locusID+")were nearly joined into the same supercontig!!");
					continue;
				}

				SuperContig tempSuper=contigB.superContig;
				Contig tempContig=tempSuper.firstContig;
				while(tempContig!=null){
					tempContig.superContig=contigA.superContig;
					tempContig=tempContig.next;
				}
				contigA.superContig.lastContig.next=tempSuper.firstContig;
				contigA.superContig.lastContig=tempSuper.lastContig;
				tempSuper.firstContig=null;
				tempSuper.lastContig=null;
			}
		}

		//NOW ASSIGN THE UNPAIRED CONTIGS TO SUPERCONTIGS (WITH ONE MEMBER EACH)
		for(int i=0; i<nContigs; i++){
			System.out.print("\rCounting reads in contig "+i);
			if(contigs[i].superContig==null){
				if(i>=contigs.length || nSuperContigs>=superContigs.length){
					//System.out.println(i+"\t"+contigs.length+"\t"+nSuperContigs+"\t"+superContigs.length);
				}
				contigs[i].superContig=superContigs[nSuperContigs];
				//NOTE THAT LOCUSID IS ASSIGNED BELOW WHEN READS IN EACH SUPERCONTIG ARE COUNTED...
				nSuperContigs++;
			}
		}

		//COUNT THE NUMBER OF READS BELONGING TO EACH CONTIG AND SUPERCONTIG
		BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(readFileStem+".ass") ) ));
		String tempS=br.readLine();
		int contigID=-1;
		while(tempS!=null){
			contigID=Integer.parseInt(tempS.split("_con")[1].split("_")[0]);
			
			//contigs[contigID-1].nReads++;
			if(contigs[contigID-1].superContig==null || superContigs[contigs[contigID-1].superContig.superContigID-1]==null){
				System.out.println(contigID+"\t");
				System.out.println(contigs[contigID-1].superContig+"\t");
				System.out.println(contigs[contigID-1].superContig.superContigID+"\t");
				System.out.println(superContigs[contigs[contigID-1].superContig.superContigID-1]);
			}
			superContigs[contigs[contigID-1].superContig.superContigID-1].nReads++;
			superContigs[contigs[contigID-1].superContig.superContigID-1].locusID=Integer.parseInt(tempS.substring(4,tempS.indexOf("_")));
			tempS=br.readLine(); //read
			tempS=br.readLine(); //qual			
			tempS=br.readLine(); //head
		}
		br.close();

		//SORT THE SUPERCONTIGS BY LOCUS AND NREADS  (NOTE: COULD BE FASTER IF LOC LOOP WAS PUT INSIDE SUPERCONTIG LOOP)
		System.out.println();
		System.out.println("\tSorting the supercontigs...");
		SuperContig sortedSCLists[] = new SuperContig[nLoci];	//each locus will have a linked list of supercontigs in descending order by numb unique Reads mapped
		for(int loc=1; loc<=nLoci; loc++){
			//create a sorted linked list of supercontigs
			SuperContig firstSC=null;
			SuperContig tempSC=null;
			for(int i=0; i<nSuperContigs; i++){
				if(superContigs[i].locusID==loc){
					if(firstSC==null){
						firstSC=superContigs[i];
					}else if(superContigs[i].nReads>firstSC.nReads){ //put in top of list
						superContigs[i].next=firstSC;
						firstSC=superContigs[i];
					}else{
						tempSC=firstSC;
						while(tempSC.next!=null){
							if(superContigs[i].nReads>tempSC.next.nReads){
								superContigs[i].next=tempSC.next;
								tempSC.next=superContigs[i];
								break;
							}
							tempSC=tempSC.next;
						}
						if(tempSC.next==null){
							tempSC.next=superContigs[i];
						}
					}
				}
			}
			sortedSCLists[loc-1]=firstSC;	//connect the list to the master array
			
		}



		//NOW OPEN THE OUT READS FILE AND REFINE THE ASSEMBLY
		
		//Store reads in linked list, sorted by alignPos and find out min and max alignPos
		//totalReads=0;	//across all loci
		//int uniqueReads=0;	//across all loci

		bw = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(readFileStem+"_assembly.fasta") ) ));	//assembly
		BufferedWriter bwConSeqs = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(readFileStem+"_conSeqs.fasta") ) ));
		BufferedWriter bwPoly = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(readFileStem+"_polys.txt") ) ));	//outfile for testing assmembly output

		//BufferedWriter bwACE = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(readFileStem+".ace") ) ));	//outfile for testing assembly output
		for(int loc=1; loc<=nLoci; loc++){
			br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(readFileStem+".ass") ) ));
			tempS=br.readLine();
			String seq;
			String locStart = ">loc"+loc+"_";
			int currReadIDX=-1;
			int currAlignPos=0;
			Contig currCon;
			SuperContig currSup;
			AlignedRead topRead=null;
			AlignedRead tempRead=null;
			Integer dupCount;
			boolean skip=false;
			int nDisabled=0;
			int nReadsMatchingLocus=0;
			int nSuperContigsThisLocus=0;

			int tempReadID;

			int conID;
			int supID;			
			//uniqueReads=0;
			//totalReads=0;

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////// A S S E M B L Y   L O A D I N G ////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// System.out.print("\rLocus "+loc+": loading mapped reads...                          ");
			//LOAD IN THE READS AND STORE IN ALIGNEDREAD LINKED LIST, SORTING AS YOU GO

			while(tempS!=null){
				skip=false;
				if(!tempS.startsWith(locStart)){		//skip reads that don't match to this locus....
					skip=true;
				}
				if(skip==false){
					nReadsMatchingLocus++;
				}
				tempReadID=Integer.parseInt(tempS.split("_read")[1].split("_")[0]);
				if(!skip && kmerRefs[tempReadID-1].disabled){		//skip reads that come from disabled kmers....
					skip=true;
					nDisabled++;
				}				
				
				if(skip){
					tempS=br.readLine();	//read
					tempS=br.readLine();	//qual
					tempS=br.readLine();	//head
					continue;					
				}

				currAlignPos=Integer.parseInt(tempS.split("_aPos")[1].split("_")[0]);
				currCon=contigs[Integer.parseInt(tempS.split("_con")[1].split("_")[0])-1];
				currSup=currCon.superContig;
				currReadIDX=Integer.parseInt(tempS.split("_read")[1].split("_")[0]);
				AlignedRead newRead = new AlignedRead();
			
				newRead.seq=br.readLine();
				newRead.qual=br.readLine();

				//dupCount = (Integer)dupMap.get(loc+"_"+currCon.contigID+"_"+currAlignPos+"_"+newRead.seq.length());

				//totalReads+=dupCount;
				//uniqueReads+=1;

				newRead.head=tempS+"dup"+dupCounts[currReadIDX-1]+"_sup"+currSup.superContigID;
				newRead.alignPos=currAlignPos;

				conID=currCon.contigID;
				supID=currSup.superContigID;
				topRead=superContigs[supID-1].firstRead;
				if(topRead==null){									//not in list yet
					topRead=newRead;
					superContigs[supID-1].firstRead=topRead;
				}else if(topRead.alignPos>newRead.alignPos){		//goes in top of list
					newRead.next=topRead;
					topRead.prev=newRead;
					topRead=newRead;
					superContigs[supID-1].firstRead=topRead;
				}else{												//goes in middle
					tempRead = topRead;
					while(tempRead.next!=null){
						if(tempRead.next.alignPos>newRead.alignPos || (tempRead.next.alignPos==newRead.alignPos && (tempRead.next.alignPos+tempRead.next.seq.length() >= newRead.alignPos+newRead.seq.length())) ){	//goes in middle of list
							newRead.next=tempRead.next;
							tempRead.next.prev=newRead;
							tempRead.next=newRead;
							newRead.prev=tempRead;
							break;
						}
						tempRead=tempRead.next;
					}
					if(tempRead.next==null){
						tempRead.next=newRead;						//goes at end of line
						newRead.prev=tempRead;
					}
				}
				tempS=br.readLine();	//head
			}
			br.close();	

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////// F I N A L I Z I N G  A S S E M B L E S /////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			AlignedRead currRead=null;
			AlignedRead refRead=null;
			int newAlignPos=-99999;
			int readLen;

			int minPos;
			int maxPos;
			int currLen;
			int newMinPos;
			int newMaxPos;
			int newPosScores[];
			int index;
			String kmer;
			int bestNewPos;
			int bestNewPosScore;

			int quickPosA;
			int quickPosB;
			int quickDif;
			boolean quickPosWorked;

			int index1=-1;
			int index2=-1;
			int breakpoint=-1;
			AlignedRead matchedRead1=null;
			AlignedRead matchedRead2=null;

			currSup=sortedSCLists[loc-1];		
			while(currSup!=null){
				if(currSup.firstRead==null || currSup.firstRead.next==null || currSup.firstRead.next.next==null){	//NOTE THAT AT LEAST THREE READS ARE REQUIRED FOR A SUPERCONTIG TO BE UTILIZED
					currSup=currSup.next;
					continue;
				}

				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////// A D J U S T I N G  A S S E M B L I E S /////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				System.out.print("\rLocus "+loc+": refining assemblies...                          ");

				//ADJUST READ START POS, AND TRIM OR INSERT GAPS IF NECESSARY

				//(trim or insert gaps on left side first.)
				currRead=currSup.firstRead.next;

				while(currRead!=null){
					currLen=currRead.seq.length();										
					//adjust alignment position if necessary

					//first attempt to see if good match to previous read
					quickPosWorked=false;
					if(currLen>30){
						quickPosA=currRead.prev.seq.indexOf(currRead.seq.substring(0,15));
						quickPosB=currRead.prev.seq.indexOf(currRead.seq.substring(currLen/2,currLen/2+15));
						if(currLen/2==quickPosB-quickPosA){
							currRead.alignPos=currRead.prev.alignPos+quickPosA;
							quickPosWorked=true;
						}
					}

					if(!quickPosWorked){	//going to have to do it the slow way...
						//what is maximum number of matches you can get to last 5 reads
						//find min and max pos of last 5 references
						minPos=99999;
						maxPos=-99999;
						refRead=currRead.prev;
						for(int i=0; i<5; i++){
							if(refRead.alignPos<minPos){
								minPos=refRead.alignPos;
							}
							if(refRead.alignPos+refRead.seq.length()>maxPos){
								maxPos=refRead.alignPos+refRead.seq.length();
							}
							refRead=refRead.prev;
							if(refRead==null){break;}
						}
					

						newMinPos=minPos-currLen;
						newMaxPos=maxPos;
						newPosScores = new int[newMaxPos-newMinPos];	//indexes are relative to newMinPos

						for(int beg=0; beg<currLen; beg++){
							for(int end=beg+alignK; end<currLen; end++){
								kmer=currRead.seq.substring(beg,end);
								refRead=currRead.prev;

								for(int i=0; i<5; i++){
									index=refRead.seq.indexOf(kmer);
									if(index>=0){
										newPosScores[refRead.alignPos+index-newMinPos-beg]+=end-beg;
									}
									refRead=refRead.prev;
									if(refRead==null){break;}
								}
							}
						}
					
						//determine best choice for offset
						bestNewPos=0;
						bestNewPosScore=0;
						for(int i=0; i<newPosScores.length; i++){
							if(newPosScores[i]>bestNewPosScore){
								bestNewPosScore=newPosScores[i];
								bestNewPos=i+newMinPos;
							}
						}
					
						//correct bestOffset to reflect relative to currPos
						newAlignPos=-99999;
						if(bestNewPosScore>0){
							newAlignPos=bestNewPos;
						}
						currRead.alignPos=newAlignPos;
					}

					if(newAlignPos<-9999){
						//remove from alignment...
						currRead.removed=true;						
						currRead.prev.next=currRead.next;
						if(currRead.next!=null){currRead.next.prev=currRead.prev;}
						currSup.nReads--;
						//currRead.head+="_NEED_TO_REMOVE_FROM_ALIGNMENT?";
					}else{
						//look for bad end ..left first
						index1=-1;
						index2=-1;
						breakpoint=-1;
						matchedRead1=null;
						matchedRead2=null;

						//FIRST IDENTIFY THE BREAKPOINT (THE FIRST POSITION IN THE CURRREAD THAT MATCHES EXACTLY TO ONE OF THE 5 PREVIOUS READS)
						for(int i=0; i<=currLen-alignK; i++){
							kmer=currRead.seq.substring(i,i+alignK);
							refRead=currRead.prev;
							for(int j=0; j<5; j++){
								index=refRead.seq.indexOf(kmer);
								//if(refRead.seq.substring(currRead.alignPos-refRead.alignPos).equals(kmer)){
								//	index=currRead.alignPos-refRead.alignPos;
								//}
								if(index>=0){// && !refRead.modified){			//match found
									//System.out.println("\tMATCH FOUND for i="+i+" j="+j+" index="+index+"\t"+(currRead.alignPos+i)+":"+(refRead.alignPos+index)+"\t"+currRead.alignPos+"\t"+refRead.alignPos);
									if(i==0 && index1==-1){			//leftmost kmer
										index1=index;	//store the results of the match with the first kmer
										matchedRead1=refRead;
									}
									if((currRead.alignPos+i)==(refRead.alignPos+index)){
										//System.out.println("\tBREAKING OUT");
										index2=index; //store the position of the match that lines up properly	(The first one)
										breakpoint=i;							
										matchedRead2=refRead;
										break;
									}
								}
								if(refRead.modified){j--;}
								refRead=refRead.prev;
								if(refRead==null){break;}
							}
							if(index2>=0){
								break;
							}
						}

						if(index1>=0 && index1==index2){
							//System.out.println("\t//NO CHANGE NEEDED");
						}else if(index1>=0 && index1<index2){
							//System.out.println("\t//SHIFT NEEDED");
							int gapLength=currRead.alignPos-(matchedRead1.alignPos+index1);
							//System.out.println("\tGap Length="+gapLength);
							//....change things....
							String newSeq=currRead.seq.substring(0,breakpoint);
							String newQual=currRead.qual.substring(0,breakpoint);
							for(int i=0; i<gapLength; i++){
								newSeq+="-";
								newQual+="!";
							}
							newSeq+=currRead.seq.substring(breakpoint);
							newQual+=currRead.qual.substring(breakpoint);
							
							currRead.seq=newSeq;
							currRead.qual=newQual;
							currRead.alignPos=currRead.alignPos-gapLength;
							currRead.modified=true;
						}else if(index1==-1 && index2>=0){
							//System.out.println("\t//TRIM NEEDED IF MISMATCHES EXCEEDS 15%");
							int nMismatches=0;
							int nCompared=0;
							for(int i=1; i<=breakpoint; i++){
								nCompared++;
								//System.out.println("\t\t"+breakpoint+"..."+i+"..."+index2);
								if(i>index2){break;}
								if(currRead.seq.charAt(breakpoint-i)!=matchedRead2.seq.charAt(index2-i)){
									nMismatches++;
								}
							}
							//System.out.println("\tnMismatches="+nMismatches+"\tnCompared="+nCompared);
							if(nMismatches>1 && nMismatches/(double)nCompared>0.15){
								//System.out.println("\tYep, trim needed...");
								///....INDICATE TRIM VIA LOWERCASE....
								//System.out.println(currRead.seq+" "+breakpoint+" "+currRead.alignPos);
								currRead.seq=currRead.seq.substring(0,breakpoint).toLowerCase()+currRead.seq.substring(breakpoint);
								//NO REASON TO CHANGE SINCE MASKED currRead.alignPos=currRead.alignPos+breakpoint;
								//System.out.println(currRead.seq+" "+breakpoint+" "+currRead.alignPos);
								currRead.modified=true;
							}
						}else if(index1==-1 && index2==-1 && currRead.prev!=null){
							//System.out.println("\t//NEED TO THROW AWAY READ, DOES NOT MATCH ANYWHERE!!");
						}else{
							//System.out.println("\t//SHOULD NOT BE HERE LEFT");
							//System.out.println(loc+"\t"+currSup.superContigID+"\t"+currRead.head.split("_read")[1].split("_")[0]+"\t"+index1+"\t"+index2+"\t"+breakpoint+"\t"+currRead.prev+"\t"+matchedRead1+"\t"+matchedRead2);
						}
					}

					currRead=currRead.next;
				}
	
				//FIX BAD RIGHT ENDS...
				//find the end of the list
				currRead=currSup.firstRead;
				while(currRead.next!=null){
					currRead=currRead.next;
				}
				int posOfEndOfLastRead=currRead.alignPos+currRead.seq.length();
				currRead=currRead.prev;	//dont mess with the last read because you won't have anything to compare it to.
				while(currRead!=null){
					currLen=currRead.seq.length();
					if(currRead.alignPos+currLen>posOfEndOfLastRead){	//this prevents end reads from being trimmed if they stick out beyond the end of the last read in the supercontig
						currRead=currRead.prev;
						continue;
					}

					if(currRead.alignPos==-99999){
						currRead.removed=true;	
						currSup.nReads--;
						if(currRead.prev!=null){currRead.prev.next=currRead.next;}
						currRead.next.prev=currRead.prev;
					}else{
						//look for bad end ..right first
						index1=-1;
						index2=-1;
						breakpoint=-1;
						matchedRead1=null;
						matchedRead2=null;

						//FIRST IDENTIFY THE BREAKPOINT (THE FIRST POSITION IN THE CURRREAD THAT MATCHES EXACTLY TO ONE OF THE 25 PREVIOUS READS)
						for(int i=currLen-alignK; i>0; i--){
							kmer=currRead.seq.substring(i,i+alignK);
							refRead=currRead.next;
							for(int j=0; j<25; j++){
								if(refRead.alignPos+refRead.seq.length()<currRead.alignPos+currLen){	//end of ref must be at or beyond end of currRead else skip
									if(refRead.modified){j--;}
									refRead=refRead.next;
									if(refRead==null){break;}
									continue;									
								}
								index=refRead.seq.indexOf(kmer);
								if(index>=0){// && !refRead.modified){			//match found
									//System.out.println("\tMATCH FOUND for i="+i+" j="+j+" index="+index+"\t"+(currRead.alignPos+i)+":"+(refRead.alignPos+index)+"\t"+currRead.alignPos+"\t"+refRead.alignPos);
									if(i==currLen-alignK && index1==-1){			//rightmost kmer
										index1=index;	//store the results of the match with the first kmer
										matchedRead1=refRead;
									}
									if((currRead.alignPos+i)==(refRead.alignPos+index)){
										//System.out.println("\tBREAKING OUT");
										index2=index; 			//store the position of the match that lines up properly	(The first one)
										breakpoint=i+alignK;	//first position mismatch occurs							
										matchedRead2=refRead;
										break;
									}
								}
								if(refRead.modified){j--;}
								refRead=refRead.next;
								if(refRead==null){break;}
							}
							if(index2>=0){
								break;
							}
						}
						//System.out.println("\t"+currRead.head.split("_read")[1].split("_")[0]+" index1:"+index1+" index2:"+index2+" breakpoint:"+breakpoint);
						if(index1>=0 && index1==index2){
							//System.out.println("\t//NO CHANGE NEEDED");
						}else if(breakpoint==-1){
							//could not find breakpoint, removing sequence
							currRead.removed=true;
							if(currRead.prev!=null){currRead.prev.next=currRead.next;}
							currRead.next.prev=currRead.prev;
							
							currSup.nReads--;

						}else if(index1>=0 && index2<index1){
							int gapLength=matchedRead1.alignPos+index1 - (currRead.alignPos+currLen-alignK);

							//System.out.println("\t//SHIFT NEEDED");														
							//System.out.println("\tGap Length="+gapLength);
							//....change things....
							//System.out.println(currRead.seq.length()+"\t"+breakpoint+"\t"+currRead.seq);
							String newSeq=currRead.seq.substring(0,breakpoint);
							String newQual=currRead.qual.substring(0,breakpoint);
							if(gapLength>=0){
								for(int i=0; i<gapLength; i++){
									newSeq+="-";
									newQual+="!";
								}	
								newSeq+=currRead.seq.substring(breakpoint);
								newQual+=currRead.qual.substring(breakpoint);
								
							}
							currRead.seq=newSeq;
							currRead.qual=newQual;
							
							//NO NEED TO CHANGE ALIGN POS currRead.alignPos=currRead.alignPos-gapLength;
							currRead.modified=true;
						}else if(index1==-1 && index2>=0){
							//System.out.println("\t//TRIM NEEDED IF MISMATCHES EXCEEDS 15%");
							int nMismatches=0;
							int nCompared=0;
							for(int i=0; i<currRead.seq.length()-breakpoint; i++){
								nCompared++;
								//System.out.println("\t\t"+breakpoint+"..."+i+"..."+index2);
								if(index2+alignK+i>=matchedRead2.seq.length()-1){break;}
								if(currRead.seq.charAt(breakpoint+i)!=matchedRead2.seq.charAt(index2+alignK+i)){
									nMismatches++;
								}
							}
							//System.out.println("\tnMismatches="+nMismatches+"\tnCompared="+nCompared);
							if(nMismatches>1 && (nCompared==0 || nMismatches/(double)nCompared>0.15)){
								//System.out.println("\tYep, trim needed...");
								///....INDICATE TRIM VIA LOWERCASE....
								//System.out.println(currRead.seq+" "+breakpoint+" "+currRead.alignPos);
								currRead.seq=currRead.seq.substring(0,breakpoint)+currRead.seq.substring(breakpoint).toLowerCase();
								//NO NEED TO CHANGE currRead.alignPos
								//System.out.println(currRead.seq+" "+breakpoint+" "+currRead.alignPos);
								currRead.modified=true;
							}
						}else if(index1==-1 && index2==-1 && currRead.prev!=null){
							//System.out.println("\t//NEED TO THROW AWAY READ, DOES NOT MATCH ANYWHERE!!");
						}else{
							//System.out.println("\t//SHOULD NOT BE HERE RIGHT");
							//System.out.println(loc+"\t"+currSup.superContigID+"\t"+currRead.head.split("_read")[1].split("_")[0]+"\t"+index1+"\t"+index2+"\t"+breakpoint+"\t"+currRead.prev+"\t"+matchedRead1+"\t"+matchedRead2);				
						}
					}

					currRead=currRead.prev;
				}			
				//remove reads marked as removed
				currRead=currSup.firstRead;
				while(currRead!=null){
					if(currRead.removed){
						//System.out.println("REMOVED READ DETECTED!!");
						//System.out.println("\t"+currRead.removed+"\t"+currRead.prev+"\t"+currRead.next+"\t"+currSup.firstRead);
						currSup.nReads--;
						if(currRead==currSup.firstRead){			//remove top node
							currSup.firstRead=currRead.next;
							if(currRead.next!=null){currRead.next.prev=null;}
						}else if(currRead.next==null){			//remove bottom node
							if(currRead.prev!=null){currRead.prev.next=null;}
						}else{	//in middle
							currRead.prev.next=currRead.next;
							currRead.next.prev=currRead.prev;
							currRead.prev=null;
						}
					}
					currRead=currRead.next;
				}

				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////// C O N S E N S U S  S E Q U E N C E S ///////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				 System.out.print("\rLocus "+loc+": computing consensus sequences...                          ");
				if(currSup.firstRead==null){
					currSup.conSeq="N";
				}else{
					nSuperContigsThisLocus++;
					currRead=currSup.firstRead;
					//first identify beg and end of alignment
					currSup.minAlignPos=currSup.firstRead.alignPos;
					currSup.maxAlignPos=currSup.firstRead.alignPos+currSup.firstRead.seq.length();

					//Determine min a max align position
					while(currRead!=null){
						if(currRead.removed){
							System.out.println("SHOULD NOT BE ANY REMOVED READS !!!!!@@@@");
							currRead=currRead.next;
							continue;
						}
						if(currRead.alignPos<currSup.minAlignPos){
							currSup.minAlignPos=currRead.alignPos;
						}
						if(currRead.alignPos+currRead.seq.length()>currSup.maxAlignPos){
							currSup.maxAlignPos=currRead.alignPos+currRead.seq.length();
						}
						currRead=currRead.next;
					}

					int alignLen=currSup.maxAlignPos-currSup.minAlignPos;

					int coverage[][] = new int[6][alignLen]; 		//A,T,C,G,N
					currRead=currSup.firstRead;				
					while(currRead!=null){
						if(!currRead.removed){
							for(int i=0; i<currRead.seq.length(); i++){
								switch(currRead.seq.charAt(i)){
									case 'A': coverage[0][currRead.alignPos-currSup.minAlignPos+i]++; break;
									case 'T': coverage[1][currRead.alignPos-currSup.minAlignPos+i]++; break;
									case 'C': coverage[2][currRead.alignPos-currSup.minAlignPos+i]++; break;
									case 'G': coverage[3][currRead.alignPos-currSup.minAlignPos+i]++; break;
									case '-': coverage[4][currRead.alignPos-currSup.minAlignPos+i]++; break;
									case 'N': coverage[5][currRead.alignPos-currSup.minAlignPos+i]++; break;
									//note that trimmed bases will be lower case and thus ignored
								}
							}
						}else{
							System.out.println("SHOULD NOT BE ANY REMOVED READS *&^");
						}
						currRead=currRead.next;
					}
				
					//call consensus bases and identify polymorphic sites
					currSup.polymorphic = new boolean[alignLen];
					currSup.nPolymorphic=0;
					currSup.conSeq="";
					for(int i=0; i<alignLen; i++){
						int sumCov=0;
						int maxCov=0;
						char maxBase='X';
						for(int j=0; j<5; j++){
							sumCov+=coverage[j][i];
							if(coverage[j][i]>maxCov){
								maxCov=coverage[j][i];
								switch(j){
									case 0: maxBase='A'; break;
									case 1: maxBase='T'; break;
									case 2: maxBase='C'; break;
									case 3: maxBase='G'; break;
									case 4: maxBase='-'; break;
									case 5: maxBase='N'; break;
								}
							}
						}

						if(sumCov==0){	//no sequences overlap with this base
							currSup.conSeq+="N";
						}else if(sumCov<5){	//coverage is too low to be sure that another allele doesn't exist, so cannot call base
							//> pbinom(0,5,0.5)
							//[1] 0.03125			//good evidence that second allele does not exist
							//> pbinom(0,4,0.5)
							//[1] 0.0625			//can't be sure that second allele does not exist
							currSup.conSeq+=callBase(coverage,i).toLowerCase();
						}else if(maxCov==sumCov){	//no disagreements, so just call base
							currSup.conSeq+=maxBase;						
						}else if(pBinom(sumCov-maxCov,sumCov,0.1)>=0.05){	//sequencing error could explain polymorphism, so just call most common base
							//****NOTE ADJUSTMENT MAY BE NEED FOR POLYPLOIDS!
							//> pbinom(25,70,0.1,lower.tail=FALSE)
							//[1] 1.319592e-09							//very unlikely that 25 sequencing errors seen in 70 reads
							//> pbinom(5,70,0.1,lower.tail=FALSE)
							//[1] 0.7127784								//reasonable to expect 5 sequencing errors in 70
							currSup.conSeq+=maxBase;
						}else{		//sequencing error cannot explain polymorphism, so need to consider this site as polymorphic due to multiple alleles
							currSup.polymorphic[i]=true;
							currSup.nPolymorphic++;
							currSup.conSeq+=callBase(coverage,i);
						}
					
					}	
				}
				
				//write the consensus sequence
				bwConSeqs.write(">L"+loc+"."+nSuperContigsThisLocus+"\n");
				bwConSeqs.write(currSup.conSeq+"\n");
				bwConSeqs.flush();

				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////// W R I T E  A S S E M B L Y /////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				 System.out.print("\rLocus "+loc+": writing assembly...                          ");

				tempRead=currSup.firstRead;

				while(tempRead!=null){
					if(tempRead.removed){
						System.out.println("SHOULD NOT BE ANY REMOVED READS !@#$");
						tempRead=tempRead.next;
						continue;
					}

					currAlignPos=tempRead.alignPos;
					currAlignPos=currAlignPos-currSup.minAlignPos;	//adjust the alignment position

					int oldAlignPos=Integer.parseInt(tempRead.head.split("_aPos")[1].split("_")[0]);
					int currReadID=Integer.parseInt(tempRead.head.split("_read")[1].split("_")[0]);

					tempRead.head=tempRead.head.replace("_aPos"+oldAlignPos,"_aPos"+currAlignPos);
					
					tempRead.alignPos=currAlignPos;
					
					//write fasta alignment format
					bw.write(">L"+loc+"."+nSuperContigsThisLocus+tempRead.head.substring(tempRead.head.indexOf("_"))+"\n");//+"_aLen"+(currSup.maxAlignPos-currSup.minAlignPos+1)+"\n");
					bw.write(tempRead.seq+"\n");
					bw.write(tempRead.qual+"\n");



					tempRead=tempRead.next;		
				}				
				
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////// P O L Y M O R P H I C   S I T E S //////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				 System.out.print("\rLocus "+loc+": identifying polymorphic sites...                          ");

				//first, make a linked list of all bases found in polymorphic sites, sorted in two ways 1) by readID, 2) by loc.sup.read.pos
				//PolyBase topPolyByRead=null;
				//PolyBase topPolyByLocus=null;
				//PolyBase bottomPolyByLocus=null;
				//PolyBase tempPolyBase=null;
				
				//obtain the polymorphic site positions
				int nPoly=currSup.nPolymorphic;
				if(nPoly>0){

					int polyPos[] = new int[nPoly];
					int counter=0;
					for(int i=0; i<currSup.polymorphic.length; i++){
						if(currSup.polymorphic[i]){
							polyPos[counter]=i;
							counter++;
						}
					}
					//System.out.println("\t"+counter);
				
					currRead=currSup.firstRead;
					int beg=-1;
					int end=-1;

					while(currRead!=null){
						beg=currRead.alignPos;
						end=currRead.alignPos+currRead.seq.length()-1;
						for(int i=0; i<nPoly; i++){
							if(beg<=polyPos[i] && end>=polyPos[i]){
								//found overlap, write the polymorphism information to the file
								bwPoly.write(loc+"\t"+nSuperContigsThisLocus+"\t"+currRead.head.split("_read")[1].split("_")[0]+"\t"+polyPos[i]+"\t"+(polyPos[i]-beg)+"\t"+currRead.seq.charAt((polyPos[i]-beg))+"\t"+currRead.qual.charAt((polyPos[i]-beg))+"\n");
							}
						}
						currRead=currRead.next;
					}
					
				}			

				 System.out.print("\rLocus "+loc+": cleaning up...                          ");

				currSup.firstRead=null;
				//System.gc();				//this slows things down a lot...
				currSup=currSup.next;
			}
			bw.flush();
			bwConSeqs.flush();
			bwPoly.flush();
		}
		bw.flush();
		bw.close();
		bwConSeqs.flush();
		bwConSeqs.close();
		bwPoly.flush();
		bwPoly.close();
		
		 //System.out.println("\t"+((new Date().getTime()-begTime)/1000)+" seconds elapsed so far.");

		 //remove temp files
		 new File("../I"+ind+"/I"+ind+".ass").delete();
		 new File("../I"+ind+"/I"+ind+".con").delete();
		 new File("../I"+ind+"/I"+ind+".cov").delete();
		 new File("../I"+ind+"/I"+ind+".merCount").delete();
					
		// long endTime = new Date().getTime();
		// System.out.println("\nSeconds required: "+(endTime-begTime)/1000.0);
      }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>> MESSAGE:"+ioe.getMessage());}
	}

  static void hashLocalKmers(String read, Read matchedRead, int round){
		Kmer tempKmer=null;
		String kmerSeq;
		Kmer currKmer;
		
		boolean badKmer;
		for(int i=0; i<=read.length()-localK; i++){
			kmerSeq=read.substring(i,i+localK);
			badKmer=false;		
			currKmer = (Kmer)localMap.get(kmerSeq);
			if(currKmer==null){		//first time observed, create and intialize				
				if(matchedRead.contigID==-1){			//this read was matched to reference (not local) so need new contigID
					nContigs++;
					matchedRead.contigID=nContigs;
				}			

				for(int j=0; j<localK-15; j++){
					if((Integer)bad15merMap.get(kmerSeq.substring(j,j+15))!=null){
						badKmer=true;
						break;
					}
				}
				
				if(badKmer){
					continue;
				}
				
				currKmer = new Kmer();
				currKmer.parent=matchedRead.kmer;
				currKmer.locus=matchedRead.locusID;
				currKmer.ref=matchedRead.refID; 			//NEW: local kmer is assined ref based on ancesgtor kmer ref, OLD: local kmer is assigned ref based on best match found in getMatches. if this read matched by locak kmer ,then refID will be inherited from the matching kmer
				currKmer.alignPos = matchedRead.alignPos+i;
				
				currKmer.depth=currKmer.parent.depth+1;
				currKmer.coverage=1;
				currKmer.contig=matchedRead.contigID;

				currKmer.roundLastMatched=round;
				currKmer.readLastMatched=matchedRead.readID;

				nLocalKmers++;						

				//put in local hashmap
				localMap.put(kmerSeq,currKmer);

			}else{	//this local kmer was observed before, need to associate the read to a preexisting contig

				currKmer.coverage++;
				if(matchedRead.contigID==-1){				//first kmer in this read, so give it a pre-existing ID
					matchedRead.contigID=currKmer.contig;
				}else if(matchedRead.contigID==currKmer.contig){
					//do nothing, kmer already exists and kmer contig matches read contig...
				}else if(matchedRead.locusID!=currKmer.locus){
						//loci for two kmers do not match, disable read!!!
System.out.println("YOU SHOULD NOT BE HERE: read"+matchedRead.readID+" loc"+matchedRead.locusID+" loc"+currKmer.locus);
				}else{

//WILL ASSIGN CONTIG PAIRS LATER
//					if(matchedRead.contigID<currKmer.contig){	//lower number goes first
//						if(contigPairMap.get(matchedRead.contigID+"_"+currKmer.contig)==null){
//							nContigPairs++;	//first time observed
//							contigPairMap.put(matchedRead.contigID+"_"+currKmer.contig,(matchedRead.alignPos+"_"+currKmer.alignPos+"_"+i+"_"+matchedRead.locusID));
//						}
//					}else{
//						if(contigPairMap.get(currKmer.contig+"_"+matchedRead.contigID)==null){
//							nContigPairs++;	//first time observed
//							contigPairMap.put(currKmer.contig+"_"+matchedRead.contigID,(matchedRead.alignPos+"_"+currKmer.alignPos+"_"+i+"_"+matchedRead.locusID));
//						}
//					}
				}
			}

		}
  }

  static void setupNmerMap(String mer){
	if(mer.length()==N){
		Nmer nmer=new Nmer();
		nmer.kmerID=-1;	//head pointer
		nmerMap.put(mer,nmer);
		return;
	}
	setupNmerMap(mer+"A");
	setupNmerMap(mer+"T");
	setupNmerMap(mer+"C");
	setupNmerMap(mer+"G");
  }

  static void loadRefSeqs(String file){
   try{	
		BufferedReader brR = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(file) ) ));
		for(int i=0; i<nRefs; i++){
//go ahead and load all of the reference sequences in...
			refStems[i]=brR.readLine();
			BufferedReader brSS = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(refStems[i]+".seqs") ) ));
			for(int j=0; j<nLoci; j++){					//no missing loci allowed!!
				//brSS.readLine();	//skip header
				refSeqs[i][j]=brSS.readLine();

				//System.out.println(refSeqs[i][j]);
				//if(refSeqs)
				refSeqLens[i][j]=refSeqs[i][j].length();
			}
			brSS.close();
		}	
		brR.close();
   }catch(IOException ioe){System.out.println("<<!!ERROR loadSeqs()!!>> "+ioe.getMessage());}		
  }

  static void loadKmers(String filestems[]){
   try{
		for(int ref=0; ref<nRefs; ref++){
			BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(filestems[ref]+".kmers") ) ));
			String tempS=br.readLine(); //e.g. loc1based	alignPos0based	0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19	G A G T G C A G C G A G A T G C C C C T
			String chunks[];
			String chunks2[];
			Kmer currKmer=null;

			String unspacedSeq="";

			Integer mapResult;

			while(tempS!=null){
				while(tempS!=null && tempS.startsWith("///////")){	//skip dividers
					tempS=br.readLine();
				}
				if(tempS==null){break;}

				//check to make sure this particular kmer has not already been used (position and locus sensitive)
				mapResult=(Integer)spacedMap.get(tempS);			
				if(mapResult!=null){
					tempS=br.readLine();
					continue;
				}
				
				spacedMap.put(tempS,0);

				chunks=tempS.split("\t");
				currKmer=new Kmer();			//new spaced kmer
				currKmer.ref=ref+1;

				//get the locus that this kmer belongs to, e.g. 402 ... 0 based now!
				currKmer.locus=Integer.parseInt(chunks[0]);

if(skipLocus[currKmer.locus-1]){//this will be true if a file named LociToSKip.txt exists and contains this locus number
	tempS=br.readLine();		
	continue;		
}

				//get position of kmer in alignment, e.g. 486
				currKmer.alignPos = Integer.parseInt(chunks[1]);

				//get positions relative to start of this kmer, e.g. 0 1 2 4 24
				chunks2=chunks[2].split(" ");			
				currKmer.sites=new int[chunks2.length];
				for(int i=0; i<chunks2.length; i++){
					currKmer.sites[i]=Integer.parseInt(chunks2[i]);	//e.g. 24
				}

				//get characters of the kmer, e.g. A T G C
				chunks2=chunks[3].split(" ");			
				currKmer.chars=new char[chunks2.length];
				for(int i=0; i<chunks2.length; i++){
					currKmer.chars[i]=chunks2[i].charAt(0);	//e.g. A
				}

				//now connect it to the linked list
				if(spacedKmers==null){
					spacedKmers=currKmer;	//first one!
				}else{
					Kmer kmer=spacedKmers;
					while(kmer.next!=null){
						kmer=kmer.next;
					}
					kmer.next=currKmer;		//put at end of linked list
				}
				nSpacedKmers++;

				//if kmer not actually spaced, also add unspaced kmer to refMap hash map
				if(currKmer.sites[19]==19){
					unspacedSeq=chunks[3].replaceAll(" ","");
					Kmer currRefKmer = (Kmer)refMap.get(unspacedSeq);
					if(currRefKmer==null){		//first time observed, create and intialize
						currKmer.kmerID=nUnspacedSpacedKmers;

						currRefKmer = new Kmer();
						currRefKmer.ref=ref+1;
						
						refMap.put(unspacedSeq,currRefKmer);

						currRefKmer.locus=currKmer.locus;

currRefKmer.disabled=currKmer.disabled;
						currRefKmer.alignPos=currKmer.alignPos;

						//now link it to the linked list
						if(refKmers==null){				//first one in list
							refKmers=currRefKmer;
						}else{							//not first one, insert at beginning
							currRefKmer.next=refKmers;
							refKmers=currRefKmer;
						}

						//put nmers from this unspaced kmer into linked list attahced to nmer in nmerMap
						for(int x=0; x<=20-N; x++){
							Nmer nmer=(Nmer)nmerMap.get(unspacedSeq.substring(x,x+N));		//will be found since all possible nmers are in map with head nodes
							if(nmer==null){continue;}
							if(nmer.next==null || nmer.next.kmerID!=nUnspacedSpacedKmers){	//only one node per kmer for each nmer
								Nmer newNmer = new Nmer();									//create a new nmer to put in the linked list
								newNmer.kmerID=nUnspacedSpacedKmers;
								newNmer.next=nmer.next;										//put new kmer at the front of the line (just after the head node)
								nmer.next=newNmer;		
							}
						}						
						nUnspacedSpacedKmers++;
					}else{	//been observed before, disable if not in same alignment position
						if(currRefKmer.locus!=currKmer.locus || currRefKmer.alignPos!=currKmer.alignPos){
							//System.out.println("Note, Kmer matching to multiple loci / alignPos found: "+unspacedSeq);
							currRefKmer.disabled=true;
							currRefKmer.locus=-1;
							currRefKmer.alignPos=-1;						
						}
					}

				}else{
					nSpacedSpacedKmers++;
				}
				tempS=br.readLine();				
			}
			br.close();
		}

   }catch(IOException ioe){System.out.println("<<!!ERROR loadKmers()!!>> "+ioe.getMessage());}	
  }

  static void loadPoss(String filestems[]){
	try{
		nSites=0;
		for(int ref=0; ref<nRefs; ref++){
			//GET THE DIMENSIONS FROM THE FILE
			BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(filestems[ref]+".poss") ) ));
			String tempS;	
			int currSites;
			for(int i=0; i<nLoci; i++){
//go ahead and load them all in
				tempS=br.readLine();					//one line per locus, 0 indicates missing data
				currSites=tempS.split(" ").length;
				if(currSites>nSites){nSites=currSites;}	//get the maximum number of sites in any alignment
			}
			br.close();
		}
		//INITIALIZE ARRAYS
		siteToPos = new int[nLoci][nRefs][nSites];
		for(int loc=0; loc<nLoci; loc++){ 
			for(int ref=0; ref<nRefs; ref++){
				for(int site=0; site<nSites; site++){
					siteToPos[loc][ref][site]=-1;
				}
			}
		}
		
		//READ IN THE VALUES
		for(int ref=0; ref<nRefs; ref++){
			BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(filestems[ref]+".poss") ) ));
			String tempS;
			String values[];
			for(int i=0; i<nLoci; i++){
				tempS=br.readLine();					//one line per locus, 0 indicates missing data
				values=tempS.split(" ");
				if(values.length<=1){
					continue;
				}
				for(int site=0; site<values.length; site++){
					siteToPos[i][ref][site]=Integer.parseInt(values[site])-1; //subtracting one because file is 1-based (0=>NA or missing). Array needs to be 0-based (-1=>NA or missing)
				}
			}
			br.close();
		}
	}catch(IOException ioe){System.out.println("<<!!ERROR loadSitePos()!!>> "+ioe.getMessage());}	
  }

  static Read getKmer(String read, HashMap kmerMap, int K){
///TRY TO SLIDE BY SOMETHING OTHER THAN 1...COULD SPEED UP CODE SUBSTANTIALLLY
	String tempSeq;
	Kmer tempKmer;
	for(int i=0; i<=read.length()-K; i++){
		tempSeq=read.substring(i,i+K);
		tempKmer=(Kmer)kmerMap.get(tempSeq);
		if(tempKmer!=null && !tempKmer.disabled && tempKmer.coverage>=kmerCoverageThreshold){			//don't match if descendent of troublemaker
			Read matchedRead=new Read();
			matchedRead.kmer=tempKmer;
			matchedRead.alignPos=tempKmer.alignPos-i;
			matchedRead.locusID=tempKmer.locus;				//1-based
			matchedRead.readID=READID;
			matchedRead.refID=tempKmer.ref;					//this allows refID to be inherited during local kmer match...  OLD: note that for ref kmer matches, .refID will be overwritten by getMatches
			matchedRead.contigID=tempKmer.contig;			//>=0 if matched to local kmer, -1 otherwise
			matchedRead.kmerPos=i;
			return(matchedRead);
		}
	}
	return null;
  }

  static Read getKmer(String read, int shortcutArray[]){
  	Kmer kmer=spacedKmers;
	kmerPosInRead=-1;
	int matches=0;
	int offset=0;
	int readLength=read.length();
	
	char readArray[]=read.toCharArray();
	
	char char1;
	char char2;
	char char3;
	char char4;
	char char5;
	char char6;
	char char7;
	char char8;
	char char9;
	char char10;
	char char11;
	char char12;
	char char13;
	char char14;
	char char15;
	char char16;
	char char17;
	char char18;
	char char19;
	char char20;
	
	int site1;
	int site2;
	int site3;
	int site4;
	int site5;
	int site6;
	int site7;
	int site8;
	int site9;
	int site10;
	int site11;
	int site12;
	int site13;
	int site14;
	int site15;
	int site16;
	int site17;
	int site18;
	int site19;
	int site20;

	
  	while(kmer!=null){
		if(kmer.disabled || kmer.notUseful){kmer=kmer.next; continue;} //don't attempt match if troublemaker or not useful		
		if(kmer.kmerID>=0 && shortcutArray[kmer.kmerID]==0){kmer=kmer.next; continue;}	//don't attempt ruled out via nmer		

		char1=kmer.chars[0];
		char2=kmer.chars[1];
		char3=kmer.chars[2];
		char4=kmer.chars[3];
		char5=kmer.chars[4];
		char6=kmer.chars[5];
		char7=kmer.chars[6];
		char8=kmer.chars[7];
		char9=kmer.chars[8];
		char10=kmer.chars[9];
		char11=kmer.chars[10];
		char12=kmer.chars[11];
		char13=kmer.chars[12];
		char14=kmer.chars[13];
		char15=kmer.chars[14];
		char16=kmer.chars[15];
		char17=kmer.chars[16];
		char18=kmer.chars[17];
		char19=kmer.chars[18];
		char20=kmer.chars[19];
	
		site1=kmer.sites[0];
		site2=kmer.sites[1];
		site3=kmer.sites[2];
		site4=kmer.sites[3];
		site5=kmer.sites[4];
		site6=kmer.sites[5];
		site7=kmer.sites[6];
		site8=kmer.sites[7];
		site9=kmer.sites[8];
		site10=kmer.sites[9];
		site11=kmer.sites[10];
		site12=kmer.sites[11];
		site13=kmer.sites[12];
		site14=kmer.sites[13];
		site15=kmer.sites[14];
		site16=kmer.sites[15];
		site17=kmer.sites[16];
		site18=kmer.sites[17];
		site19=kmer.sites[18];
		site20=kmer.sites[19];


  		//check to see if this kmer matches anywhere within the read
 		offset=0;
		while((site20+offset)<readLength){
			matches=0;

			if(char1==readArray[site1+offset]){matches++;} 					//if(1-matches>3){offset++; continue;}
			if(char2==readArray[site2+offset]){matches++;} 					//if(2-matches>3){offset++; continue;}
			if(char3==readArray[site3+offset]){matches++;} 					//if(3-matches>3){offset++; continue;}
			if(char4==readArray[site4+offset]){matches++;} 					//if(4-matches>3){offset++; continue;}
			if(char5==readArray[site5+offset]){matches++;} 					//if(5-matches>3){offset++; continue;}
			if(char6==readArray[site6+offset]){matches++;} 					//if(6-matches>3){offset++; continue;} //waste of time before 7
			if(char7==readArray[site7+offset]){matches++;} 					if(7-matches>3){offset++; continue;}
			if(char8==readArray[site8+offset]){matches++;} 					if(8-matches>3){offset++; continue;}
			if(char9==readArray[site9+offset]){matches++;} 					if(9-matches>3){offset++; continue;}
			if(char10==readArray[site10+offset]){matches++;} 				if(10-matches>3){offset++; continue;}
			if(char11==readArray[site11+offset]){matches++;} 				if(11-matches>3){offset++; continue;}
			if(char12==readArray[site12+offset]){matches++;} 				if(12-matches>3){offset++; continue;}
			if(char13==readArray[site13+offset]){matches++;} 				if(13-matches>3){offset++; continue;}
			if(char14==readArray[site14+offset]){matches++;} 				if(14-matches>3){offset++; continue;}
			if(char15==readArray[site15+offset]){matches++;} 				if(15-matches>3){offset++; continue;}
			if(char16==readArray[site16+offset]){matches++;} 				if(16-matches>3){offset++; continue;}
			if(char17==readArray[site17+offset]){matches++;} 				if(17-matches>3){offset++; continue;}
			if(char18==readArray[site18+offset]){matches++;} 				if(18-matches>3){offset++; continue;}
			if(char19==readArray[site19+offset]){matches++;} 				if(19-matches>3){offset++; continue;}
			if(char20==readArray[site20+offset]){matches++;} 				if(20-matches>3){offset++; continue;}

  			if(matches>=matchThreshold1){
				//System.out.println("FOUND FIRST PASS MATCH!");
				Read matchedRead=new Read();
				matchedRead.kmer=kmer;
				matchedRead.alignPos=kmer.alignPos-offset;
				matchedRead.locusID=kmer.locus;				//1-based
				matchedRead.readID=READID;
				matchedRead.kmerPos=offset;
				//matchedRead.spacedKmerMatchScore=matches;
				return(matchedRead);
  			}
  			offset++;
  		}
  		kmer=kmer.next;
  	}

  	return null;
  }

  static void getMatches(String readSeq, Read read){

	int nMatches=0;
	double currScore=0;
	double bestScore=0;
	int bestRef=-1;
	double currOverlap=0;
	double bestOverlap=0;
	int offset=0;
	int posInRef=-1;
	int posInRead=-1;
	int locIndex=read.locusID-1;
	int alignPos = read.kmer.alignPos;

	for(int ref=0; ref<nRefs; ref++){
		//count the number of matches at the position specified by the location the kmer matched
		nMatches=0;
		currScore=0;
		currOverlap=0;
 
		posInRead = read.kmerPos;

		posInRef = siteToPos[locIndex][ref][alignPos];
		offset=0;	//start at kmer starting position

		while(posInRead+offset>=0 && posInRef+offset>=0){
			currOverlap++;
			if( readSeq.charAt(posInRead+offset)==refSeqs[ref][locIndex].charAt(posInRef+offset)){
				nMatches++;
			}
			offset--;
		}

		posInRead = read.kmerPos;
		posInRef = siteToPos[locIndex][ref][alignPos];
		offset=1;	//start at one position to right of kmer starting position

		while(posInRead+offset<readSeq.length() && posInRef+offset<refSeqs[ref][locIndex].length()){
			currOverlap++;
			if( readSeq.charAt(posInRead+offset)==refSeqs[ref][locIndex].charAt(posInRef+offset)){
				nMatches++;
			}
			offset++;
		}

		if(currOverlap>=100 && nMatches/currOverlap > matchThreshold2){
			currScore=nMatches/currOverlap;
			if(currScore>read.matchScore){
				//keep the previous best
				//read.matchScore2nd=read.matchScore;
				//read.overlap2nd=read.overlap;
				//read.refID2nd=read.refID;
				
				//now store the best
				read.matchScore=currScore;
				//read.overlap=currOverlap;
				read.refID=ref+1;
				//read.locusID=locIndex+1; already determined based on preliminary match
			}
		}
	}
  }

  static int adjustPos2(String seq1, String seq2, int pos1, int pos2){	//assumes pos1<pos2
	int kmerPos=-1;
	for(int i=0; i<=seq2.length()-localK; i+=1){
		kmerPos=seq1.indexOf(seq2.substring(i,i+alignK));
		if(kmerPos>=0){
			return pos1+kmerPos-i;
		}
	}
	//System.out.println("Note: could not adjust pos2!");
	return -99999;
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
	
  static double pBinom(int k,int n,double prob){
		double pbinom=0;
		for(int x=k; x<=n; x++){
			//System.out.println(""+x+"\t"+choose(n,x,true) +"\t"+ Math.log(Math.pow(prob,x)) +"\t"+ Math.log(Math.pow(1-prob,n-x))+"\t"+Math.exp(choose(n,x,true) + Math.log(Math.pow(prob,x)) + Math.log(Math.pow(1-prob,n-x)))+"\t"+pbinom);			
			pbinom+=Math.exp(choose(n,x,true) + Math.log(Math.pow(prob,x)) + Math.log(Math.pow(1-prob,n-x)));	//note that working in logs avoids overflow for factorial(>171)
			//System.out.println(""+x+"\t"+choose(n,x,true) +"\t"+ Math.log(Math.pow(prob,x)) +"\t"+ Math.log(Math.pow(1-prob,n-x))+"\t"+Math.exp(choose(n,x,true) + Math.log(Math.pow(prob,x)) + Math.log(Math.pow(1-prob,n-x)))+"\t"+pbinom);
		}
		return pbinom;
  }
  
  static String callBase(int coverage[][], int i){
						//A                   T                   C                   G 
	if(      coverage[0][i]>=1 && coverage[1][i]==0 && coverage[2][i]==0 && coverage[3][i]==0){	return "A";
	}else if(coverage[0][i]==0 && coverage[1][i]>=1 && coverage[2][i]==0 && coverage[3][i]==0){	return "T";
	}else if(coverage[0][i]==0 && coverage[1][i]==0 && coverage[2][i]>=1 && coverage[3][i]==0){	return "C";
	}else if(coverage[0][i]==0 && coverage[1][i]==0 && coverage[2][i]==0 && coverage[3][i]>=1){	return "G";
	}else if(coverage[0][i]>=1 && coverage[1][i]==0 && coverage[2][i]==0 && coverage[3][i]>=1){	return "R";	
	}else if(coverage[0][i]==0 && coverage[1][i]>=1 && coverage[2][i]>=1 && coverage[3][i]==0){	return "Y";
	}else if(coverage[0][i]==0 && coverage[1][i]==0 && coverage[2][i]>=1 && coverage[3][i]>=1){	return "S";
	}else if(coverage[0][i]>=1 && coverage[1][i]>=1 && coverage[2][i]==0 && coverage[3][i]==0){	return "W";
	}else if(coverage[0][i]==0 && coverage[1][i]>=1 && coverage[2][i]==0 && coverage[3][i]>=1){	return "K";
	}else if(coverage[0][i]>=1 && coverage[1][i]==0 && coverage[2][i]>=1 && coverage[3][i]==0){	return "M";
	}else if(coverage[0][i]==0 && coverage[1][i]>=1 && coverage[2][i]>=1 && coverage[3][i]>=1){	return "B";
	}else if(coverage[0][i]>=1 && coverage[1][i]>=1 && coverage[2][i]==0 && coverage[3][i]>=1){	return "D";
	}else if(coverage[0][i]>=1 && coverage[1][i]>=1 && coverage[2][i]>=1 && coverage[3][i]==0){	return "H";
	}else if(coverage[0][i]>=1 && coverage[1][i]==0 && coverage[2][i]>=1 && coverage[3][i]>=1){	return "V";						
	}else if(coverage[0][i]>=1 && coverage[1][i]>=1 && coverage[2][i]>=1 && coverage[3][i]>=1){	return "N";
	}else{																						return "-";
	}
	
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

class Kmer{
	int kmerID=-1;				//position in NmerKmerMatches to which this kmer belongs
	int locus=-1;				//locus to which this kmer corresponds
	int ref=-1;					//reference to which this kmer belongs...OLD: ONLY RELEVANT FOR LOCAL MATCHES WITH LONG KMERS THAT ARE NOT CHECKED WITH EXTENSIVE MATCHING
	int contig=-1;				//number indicating a cluster of reads assumed to belong to the same set of assembly (could be from differnt alleles but should not be from differnt homologs)
	int alignPos=-1;			//position this kmer begins at in the reference alignment (assumed it only comes up once!)
	int roundLastMatched=1;		//this is used to allow kmers to be removed from hash table once they have been unsuccessfully checked against the entire read pile
	int readLastMatched=1;		//this is used to allow kmers to be removed from hash table once they have been unsuccessfully checked against the entire read pile
	char chars[];				//used to store spaced kmer sequence, if relevant
	int sites[];				//used to store spaced kmer positions, if relevant
	//int nChildren=0;			//number of kmers put into hash table because of this kmer
	// int nGrandchildren=0;		//number of kmers put into has table because one of it's children's kmers
	int coverage=0;				//in how many reads is this kmer found...>1 required for use of kmer for local mapping!!
	int depth=0;		
	Kmer parent=null;			//pointer to kmer
	boolean disabled=false;		//set true if descendant of troublemaker
	boolean notUseful=false;	//true if another kmer from the same locus (different reference) is providing many more matches see getKmer code for details
	Kmer next=null;				//pointer to next kmer in main linked list
}

class Nmer{						//this class is used to make nodes for the Nmer lookup table
	int kmerID=-1;				//index of spaced kmer
	Nmer next;					//pointer to next node
}

class Read{
	int readID=-1;				//index identifying the read in the input file
	int refID=-1;				//reference to which best extensive match occurred
	int locusID=-1;				//locus to which this read belongs
	int contigID=-1;			//contig to which this read belongs...will be joined to other contigs...
	int nRedundant=0;			//number of reads with identical sequence to this one
	char orientation='X';		//'F' or 'R'
	int alignPos=-999999;		//Expected position in alignment of beginning of read relative ...only a guess based on the position of the kmer in the sequence and referece.....NOTE THAT -999999 means unplaced
	int kmerPos=-999999;		//temporary variable to store the position in the read at which the kmer matched in teh prelminary match, only used in getMatches for lining up the read with the reference sequence
	int matchStatus=0;			//indicator of how the match was made 0=no match, 1=local kmer match, 2=reference kmer match, 3=reference spaced kmer match
	double matchScore=0;		//% of overlapping bases that matched
	Kmer kmer=null;
}

class AlignedRead{
	int alignPos=-999999;
	boolean modified=false;
	boolean removed=false;
	String head="";
	String seq="";
	String qual="";
	AlignedRead next;
	AlignedRead prev;
}

class Contig{
	int contigID=-1;
	int locusID=-1;
	Contig next=null;
	SuperContig superContig=null;
}

class SuperContig{
	int superContigID=-1;
	int nReads=0;
	int locusID=-1;
	String conSeq="";
	boolean polymorphic[];
	int nPolymorphic=0;
	int minAlignPos=-99999;
	int maxAlignPos=-99999;
	AlignedRead firstRead=null;
	Contig firstContig=null;
	Contig lastContig=null;
	SuperContig next;			//used for sorted linked list of supercontigs for writing reads to file
}

