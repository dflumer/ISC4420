import java.io.*;

public class OrganizeReadFiles{

  public static void main(String[] args){
    	File[] fileList = new File("../").listFiles();

    	for(int i=0; i<fileList.length; i++){
    		String filename=(""+fileList[i]);
    		filename=filename.substring(filename.lastIndexOf("/")+1);
    		if(filename.endsWith(".fastq")){
    			String ind=filename.split("_")[0];
	    		if(!new File("../"+ind).exists()){
    				new File("../"+ind).mkdir();
    			}
    			System.out.println("../"+ind+"/"+filename);
				fileList[i].renameTo(new File("../"+ind+"/"+filename));
    		}
    	}
  }

}
