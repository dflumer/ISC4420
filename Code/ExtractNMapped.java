// java ExtractNMapped 9884 10345
// for this code, the I folder needs to contain the _assembly.fasta file
// You must run assemblies before you can run this code
import java.io.*;

public class ExtractNMapped{
	public static void main(String[] args){
	  try{
		
		int begInd = Integer.parseInt(args[0]);	//e.g. 2389
		int endInd = Integer.parseInt(args[1]); //e.g. 2455
		
		for(int ind=begInd; ind<=endInd; ind++){

			String sampleName = getSampleName(ind);
			//finally, summize the assembly itself
			String filename="../"+sampleName+"/"+sampleName+"_assembly.fasta";
			if(new File(filename).exists()){

				BufferedReader br = new BufferedReader ( new InputStreamReader ( new FileInputStream (new File(filename))));
				BufferedWriter bw = new BufferedWriter ( new OutputStreamWriter ( new FileOutputStream (new File("../"+sampleName+"/"+sampleName+"_nMapped.txt"))));
				String tempS=br.readLine();
				int superID=-1;
				int locID=-1;
				int nDups=-1;
				int prevSuperID=Integer.parseInt(tempS.substring(tempS.indexOf(".")+1,tempS.indexOf("_")));
				int prevLocID=Integer.parseInt(tempS.substring(2,tempS.indexOf(".")));
				int nMappedReads=0;
				while(tempS!=null){
					superID=Integer.parseInt(tempS.substring(tempS.indexOf(".")+1,tempS.indexOf("_")));
					locID=Integer.parseInt(tempS.substring(2,tempS.indexOf(".")));
					nDups=Integer.parseInt(tempS.split("_dup")[1].split("_")[0]);

					//if(superID>100){
					//	tempS=br.readLine();
					//	tempS=br.readLine();
					//	tempS=br.readLine();	
					//	continue;
					//}

					tempS=br.readLine();//sequence

					if(superID!=prevSuperID || locID!=prevLocID){	//new consensus sequence
						System.out.print("\r"+ind);
						bw.write(""+prevLocID+"\t"+prevSuperID+"\t"+nMappedReads+"\n");
						nMappedReads=0;

					}
					prevSuperID=superID;
					prevLocID=locID;
					nMappedReads+=nDups;
					
					tempS=br.readLine();//quals
					tempS=br.readLine();//next header
				}
				bw.write(""+prevLocID+"\t"+prevSuperID+"\t"+nMappedReads+"\n");
				bw.flush();
				bw.close();
				br.close();
			}else{System.out.println("\nThere is no assembly file to make the nMapped file from.");}
		}		
		System.out.println();
	} catch(IOException ioe){System.out.println("\n\n<<!!ERROR main()!!>>"+ioe.getMessage());}
			
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
