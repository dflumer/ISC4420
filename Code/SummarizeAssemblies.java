import java.io.*;

public class SummarizeAssemblies{
	public static void main(String[] args){
	  try{
		
		int begInd = 1;
		int endInd = 100000;

		int nLoci=0;
		for(int ind=begInd; ind<endInd; ind++){
			String filename="../I"+ind+"/I"+ind+"_nMapped.txt";
			if(new File(filename).exists()){
				BufferedReader br = new BufferedReader ( new InputStreamReader ( new FileInputStream (new File(filename))));
				String tempS=br.readLine();
				while(tempS!=null){
					int locus=Integer.parseInt(tempS.split("\t")[0]);
					if(locus>nLoci){nLoci=locus;}
					tempS=br.readLine();
				}
				br.close();
			}
		}
		System.out.println(nLoci+" loci detected.");

		String outfileStem = "../Results/Assembly";			//e.g. ../Results/P0033_AssemblySummary.txt

		int begLoc=1;
		int endLoc=nLoci;


		//determine the number of individuals and get a ind ID converter so IndID and be converted to array index
		int nInds=0;
		int indIndexes[] = new int[endInd+1];
		for(int ind=begInd; ind<=endInd; ind++){
			String sampleName = getSampleName(ind);
			String filename="../"+sampleName+"/"+sampleName+"_assembly.fasta";
			if(new File(filename).exists()){ nInds++; indIndexes[ind]=nInds-1;}else{indIndexes[ind]=-1;}
		}
			
		long nSuperContigs[][]=new long[nInds][nLoci];

		if(!new File("../Results").exists()){new File("../Results").mkdir();}
		
		int maxNSuperContigs=0;
		for(int ind=begInd; ind<=endInd; ind++){
			//first determine the number of homologs per locus
			String sampleName = getSampleName(ind);

			String filename="../"+sampleName+"/"+sampleName+"_assembly.fasta";
			if(!new File(filename).exists()){continue; }
			
			System.out.print("\r"+sampleName);

			BufferedReader br = new BufferedReader ( new InputStreamReader ( new FileInputStream (new File(filename))));
			String tempS=br.readLine();
			int lineNumber=1;
			while(tempS!=null){
				if(tempS.indexOf(".")==-1){System.out.println("\n"+lineNumber+"\n");}
				//System.out.println(tempS.indexOf("."));
				//System.out.println(tempS);


				int superID=Integer.parseInt(tempS.substring(tempS.indexOf(".")+1,tempS.indexOf("_")));
				int locID=Integer.parseInt(tempS.substring(2,tempS.indexOf(".")));
				nSuperContigs[indIndexes[ind]][locID-1]=superID;
				if(superID>maxNSuperContigs){maxNSuperContigs=superID;}

				tempS=br.readLine();//sequence
				tempS=br.readLine();//quality scores
				tempS=br.readLine();
				lineNumber+=3;
				
			}
			br.close();
		}
		System.out.println();
		System.out.println("Maximum number of supercontigs="+maxNSuperContigs);
	
		if(maxNSuperContigs>100){maxNSuperContigs=100; System.out.println("Only Writing results for first 100 supercontigs.");}
	
		long nRawReadsP[] = new long[nInds];	
		long nRawReadsM[] = new long[nInds];	
		long nRawReadsU[] = new long[nInds];	

		long nRawBasesP[] = new long[nInds];	//pre-merging
		long nRawBasesM[] = new long[nInds];	//merged
		long nRawBasesU[] = new long[nInds];	//unmerged

		long conSeqLensLower[][][] = new long[maxNSuperContigs][nInds][nLoci];
		long conSeqLensUpper[][][] = new long[maxNSuperContigs][nInds][nLoci];
		long conSeqLensLowerUnAmbig[][][] = new long[maxNSuperContigs][nInds][nLoci];
		long conSeqLensUpperUnAmbig[][][] = new long[maxNSuperContigs][nInds][nLoci];

		long nMappedBases[][][] = new long[maxNSuperContigs][nInds][nLoci];
		long nMappedReads[][][] = new long[maxNSuperContigs][nInds][nLoci];
		
		long nMappedBasesAll[] = new long[nInds];

		long conSeqLensMax[][] = new long[nInds][nLoci];					//maximum number of con seq length across all homologs (includes upper and lower case ambig and unambig)
		long conSeqLensUpperMax[][] = new long[nInds][nLoci];

		double dupRateAll[]= new double[nInds]; //duplication rate over all reads
	 	double dupRateLast[]=new double[nInds]; //duplication rate over last 1000 samples
		
		for(int ind=begInd; ind<=endInd; ind++){
			int indIndex=indIndexes[ind];
			if(indIndex==-1){continue;}

			String sampleName = getSampleName(ind);
			//first determine the number of raw reads/bases
			String filename="../"+sampleName+"/"+sampleName+"_lengths.txt";
			if(new File(filename).exists()){
				BufferedReader br = new BufferedReader ( new InputStreamReader ( new FileInputStream (new File(filename))));
				System.out.println("Determining raw reads from "+filename);
				String tempS=br.readLine();
				int mergedReadLength=0;
				while(tempS!=null){
					mergedReadLength++;
					long value=Integer.parseInt(tempS);
										
					tempS=br.readLine();

					if(tempS!=null){
						nRawReadsP[indIndex]+=2*value;
						nRawReadsM[indIndex]+=value;
						nRawBasesP[indIndex]+=value;		//will multiply by total paired read length below (e.g. 300)
						nRawBasesM[indIndex]+=value*mergedReadLength;
					}else{
						nRawReadsP[indIndex]+=value;
						nRawReadsU[indIndex]+=value;	
						nRawBasesP[indIndex]+=value;											
						nRawBasesP[indIndex]*=mergedReadLength;
						nRawBasesU[indIndex]+=value*mergedReadLength;						
					}

				}
				br.close();
			}

			//second, summarize the consensus sequences
			filename="../"+sampleName+"/"+sampleName+"_conSeqs.fasta";
			
			if(new File(filename).exists()){
				BufferedReader br = new BufferedReader ( new InputStreamReader ( new FileInputStream (new File(filename))));
				System.out.println("Summarizing the consensus sequences in "+filename);

				String tempS=br.readLine();
				while(tempS!=null){
					int locID=Integer.parseInt(tempS.substring(2,tempS.indexOf(".")));
					int supID=Integer.parseInt(tempS.substring(tempS.indexOf(".")+1));
					tempS=br.readLine(); 	//sequence

					if(supID>maxNSuperContigs){
						tempS=br.readLine();
						continue;
					}

					for(int i=0; i<tempS.length(); i++){
						if(Character.isUpperCase(tempS.charAt(i))){
							conSeqLensUpper[supID-1][indIndex][locID-1]++;
						}else{
							conSeqLensLower[supID-1][indIndex][locID-1]++;
						}
						switch(tempS.charAt(i)){
							case 'A': case 'T': case 'C': case 'G': conSeqLensUpperUnAmbig[supID-1][indIndex][locID-1]++; break;
							case 'a': case 't': case 'c': case 'g': conSeqLensLowerUnAmbig[supID-1][indIndex][locID-1]++; break;
						}
					}
					if(conSeqLensUpper[supID-1][indIndex][locID-1]>conSeqLensUpperMax[indIndex][locID-1]){
						conSeqLensMax[indIndex][locID-1]=conSeqLensUpper[supID-1][indIndex][locID-1]+conSeqLensLower[supID-1][indIndex][locID-1];
						conSeqLensUpperMax[indIndex][locID-1]=conSeqLensUpper[supID-1][indIndex][locID-1];
					}
					tempS=br.readLine();//next heADER
					
				}
				br.close();				
			}

			//finally, summize the assembly itself
			filename="../"+sampleName+"/"+sampleName+"_assembly.fasta";
			if(new File(filename).exists()){
				System.out.println("Summarzing the assembly "+filename);

				BufferedReader br = new BufferedReader ( new InputStreamReader ( new FileInputStream (new File(filename))));
				String tempS=br.readLine();
				while(tempS!=null){
					int superID=Integer.parseInt(tempS.substring(tempS.indexOf(".")+1,tempS.indexOf("_")));
					int locID=Integer.parseInt(tempS.substring(2,tempS.indexOf(".")));
					int nDups=Integer.parseInt(tempS.split("_dup")[1].split("_")[0]);

					if(superID>maxNSuperContigs){
						tempS=br.readLine();
						tempS=br.readLine();
						tempS=br.readLine();	
						continue;
					}

					tempS=br.readLine();//sequence

					nMappedReads[superID-1][indIndex][locID-1]+=nDups;
					nMappedBases[superID-1][indIndex][locID-1]+=tempS.length()*nDups;
					nMappedBasesAll[indIndex]+=tempS.length()*nDups;

					tempS=br.readLine();//quals
					tempS=br.readLine();//next header
				}
				br.close();
			}

			//DETERMINE THE DUPLICATE RATE (REDUNDANT READ RATE)
			BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../"+sampleName+"/"+sampleName+"_readsInFiles.txt") ) ));
			String tempS=br.readLine();
			int nReadsM=Integer.parseInt(tempS.split("\t")[1]); //M
			int nReadsALL=nReadsM;
			tempS=br.readLine();
			nReadsALL+=Integer.parseInt(tempS.split("\t")[1]); //U1
			tempS=br.readLine();
			nReadsALL+=Integer.parseInt(tempS.split("\t")[1]); //U2	

			int dupID=0;

			//fill in the missing numbers caused by the 20mer shortcut
			int dups2[] = new int[nReadsALL];
			br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File("../"+sampleName+"/"+sampleName+".dups") ) ));

			tempS=br.readLine();
			int count=0;
			while(tempS!=null){
				dupID=Integer.parseInt(tempS);
				count++;
				//if(count>=nReadsM){break;}
				if(dupID>count){
					//unique via 20mer, need to insert it
					dups2[count-1]=count;
				}else{
					dups2[count-1]=dupID;
					tempS=br.readLine();
				}
			}
			br.close();

			int isNew[] = new int[nReadsALL];
			for(int i=0; i<nReadsALL; i++){
				if(dups2[i]==i+1){
					isNew[i]=1;
				}
			}

			//calculate the # reads new using a sliding window of 10000 reads		
			int windowTally=0;				//calculation over just last 1000 reads
			int overallTally=0;				//calculation over all samlpes
			for(int i=0; i<nReadsALL; i++){
				overallTally+=isNew[i];
				if(i>nReadsALL-1000){
					windowTally+=isNew[i];
				}
			}
			dupRateAll[indIndex]=1-overallTally/(double)nReadsALL;
			dupRateLast[indIndex]=1-windowTally/1000.0;	


		}		
		



		//compute summaries of the summaries
		System.out.println("Writing to "+outfileStem+"_Summary.csv");

		BufferedWriter bw = new BufferedWriter (new OutputStreamWriter(new FileOutputStream(new File(outfileStem+"_Summary.csv"))));
		bw.write("Individual,SampleID,Family,Genus,Epithet,Notes,TaxonSet,nRawReads,nRawBases,nLoci125,nLoci250,nLoci500,nLoci1000,AvgLocLen,%OnTarget,AvgNHomologs,DupRateAll\n");		
		
		for(int ind=begInd; ind<=endInd; ind++){
			int indIndex=indIndexes[ind];
			if(indIndex==-1){continue;}

			String sampleName = getSampleName(ind);
			File conseqs = new File("../"+sampleName+"/"+sampleName+"_conSeqs.fasta");
			int currChar = 0, totalNs = 0, totalACTGs = 0, totalOthers = 0, linenumber = 0;
			if (conseqs.isFile()){
				BufferedReader broS = new BufferedReader (new InputStreamReader(new FileInputStream(conseqs)));
				String tempQ = broS.readLine();
				
				char c;
				while(tempQ!=null){
					if (tempQ.charAt(0) != '>'){
						currChar = 0;
						while(currChar<(tempQ.length()-1)){
							switch(Character.toUpperCase(tempQ.charAt(currChar))){
								case 'A':
								case 'T':
								case 'C':
								case 'G': totalACTGs++; break;
								case 'N': totalNs++;    break;
								default : totalOthers++; break;
							}
						
							currChar++;
						}
					}
					tempQ = broS.readLine();
					linenumber++;
				}
				//System.out.println("ATCGs: "+totalACTGs);
				//System.out.println("Ns: "+totalNs);
				//System.out.println("Others: "+totalOthers);
				broS.close();
			} else System.out.println("The file specified file could not be read. Please check if it is a valid file.");
	
//////////////
			bw.write(sampleName+",");
			bw.write(",");
			bw.write(",");
			bw.write(",");
			bw.write(",");
			bw.write(",");
			bw.write(",");
			bw.write(nRawReadsP[indIndex]+",");
			bw.write(nRawBasesP[indIndex]+",");
			int nLoci125=0;
			int nLoci250=0;
			int nLoci500=0;
			int nLoci1000=0;
			
			int sumConSeqLensMax=0;
			int sumConSeqLensUpper=0;
			int sumMappedReads=0;
			int sumMappedBases=0;
			int sumNHomologs=0;
			int sumRawBases=0;

			for(int i=0; i<nLoci; i++){
				sumConSeqLensMax+=conSeqLensMax[indIndex][i];
				sumConSeqLensUpper+=conSeqLensUpperMax[indIndex][i];
				for(int s=0; s<maxNSuperContigs; s++){
					sumMappedReads+=nMappedReads[s][indIndex][i];	
					sumMappedBases+=nMappedBases[s][indIndex][i];	
					if(nMappedReads[s][indIndex][i]>10){
						sumNHomologs++;
					}
				}

				if(conSeqLensUpperMax[indIndex][i]>=1000){
					nLoci125++;
					nLoci250++;
					nLoci500++;
					nLoci1000++;
				}else if(	conSeqLensUpperMax[indIndex][i]>=500){
					nLoci125++;
					nLoci250++;
					nLoci500++;
				}else if(	conSeqLensUpperMax[indIndex][i]>=250){
					nLoci125++;
					nLoci250++;
				}else if(	conSeqLensUpperMax[indIndex][i]>=125){
					nLoci125++;
				}


			}

			bw.write(nLoci125+",");
			bw.write(nLoci250+",");
			bw.write(nLoci500+",");
			bw.write(nLoci1000+",");
			bw.write(sumConSeqLensMax/(double)nLoci+",");
			bw.write((100.0*sumMappedBases/(double)((double)nRawBasesM[indIndex]+(double)nRawBasesU[indIndex]))+",");	 //% on target
			bw.write(sumNHomologs/(double)nLoci+",");
			bw.write(dupRateAll[indIndex]+"\n");
		}
		bw.flush();
		bw.close();

	}catch(IOException ioe){System.out.println("\n\n<<!!ERROR main()!!>>"+ioe.getMessage());}	
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

}
