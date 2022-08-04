import java.io.*;

public class RecoverReadsInFiles{
  public static void main(String[] args){
      try{

      	for(int i=1; i<100000; i++){
      		String fileStem="../I"+i+"/I"+i;
      		if(new File(fileStem+"_M.fastq").exists()){
      			System.out.println("Counting lines in I"+i);
				BufferedWriter bw = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(fileStem+"_readsInFiles.txt") ) ));
      			BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(fileStem+"_M.fastq") ) ));
      			String tempS=br.readLine();
      			int nLinesM=0;
      			while(tempS!=null){
      				tempS=br.readLine();
      				nLinesM++;
      			}
      			br.close();
      			bw.write(fileStem+"_M.fastq"+"\t"+nLinesM/4+"\n");

				int nLinesU1=0;
      			br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(fileStem+"_U1.fastq") ) ));
      			tempS=br.readLine();
      			while(tempS!=null){
      				tempS=br.readLine();
      				nLinesU1++;
      			}
      			br.close();
      			bw.write(fileStem+"_U1.fastq"+"\t"+nLinesU1/4+"\n");

				int nLinesU2=0;
      			br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(fileStem+"_U2.fastq") ) ));
      			tempS=br.readLine();
      			while(tempS!=null){
      				tempS=br.readLine();
      				nLinesU2++;
      			}
      			br.close();
      			bw.write(fileStem+"_U2.fastq"+"\t"+nLinesU2/4+"\n");
      			bw.flush();
      			bw.close();

      			if(nLinesU1!=nLinesU2){
      				System.out.println("WARNING: nLinesU1!=nLinesU2 for I"+i);
      			}
      		}
      	}

      }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>>"+ioe.getMessage()+" Caution: Individual May Not Exist");}
  }
}