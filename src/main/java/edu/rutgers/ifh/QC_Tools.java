package edu.rutgers.ifh;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class QC_Tools {
	
	public String str_Cluster_User_Name = "";
	public String str_Cluster_User_Password = "";
	public String str_Cluster_Host_Name = "";
	public int int_Cluster_Host_port = 0;
	
	public String str_DB_Table = "";
	public String str_DB_Host_Name = "";
	public int int_DB_port = 0;
	
	public String str_Vcf2Sql_Location = "";
	
	public String str_job_name = "";
	public String str_FastQC_Location = "";
	public String str_Trimmomatic_Location = "";
	public String str_Picard_Location = "";
	public String str_BWA_Location = "";
	public String str_Samtools_Location = "";
	public String str_Gatk_Location = "";
	public String str_SnpEff_Location = "";
	public String str_Bedtools_Location = "";
	public String str_hg38_Location = "";
	public String str_dbsnp_Location = "";
	public String str_hapmap_Location = "";
	public String str_mills_Location = "";
	public String str_omni_Location = "";
	public String str_phase1_Location = "";
	public String str_r1_filtered_trimmed_Location = "";
	public String str_r1_unfiltered_Location = "";
	public String str_r2_filtered_trimmed_Location = "";
	public String str_r2_unfiltered_Location = "";
	public String str_TruSeq3_Location = "";

	
	public boolean mFastqc = false;
	public boolean mTrimmomatic = false;
	public boolean mPicard = false;
	public boolean mBwa = false;
	public boolean mSamtools = false;
	public boolean mGatk = false;
	public boolean mSnpEff = false;
	public boolean mBedtools = false;
 	
	InputStream inputStream;
	
	QC_Tools() {
		
	}
	
	void LoadPropertiesFiles() throws IOException {
		try {
			Properties prop = new Properties();
			String propFileName = "config.properties";
 
			inputStream = getClass().getClassLoader().getResourceAsStream(propFileName);
 
			if (inputStream != null) {
				prop.load(inputStream);
				
				str_job_name = prop.getProperty("Job.Name");
				
				str_Cluster_Host_Name = prop.getProperty("Cluster.Host.Name");
				
				str_DB_Host_Name = prop.getProperty("DB.Host.Name");
				int_DB_port = Integer.valueOf(prop.getProperty("DB.Port"));
				str_DB_Table = prop.getProperty("DB.Table.Name");
				
				str_Vcf2Sql_Location = prop.getProperty("Vcf2SQL.Path");
				 
				// get the property value and print it out
				mFastqc = Boolean.parseBoolean(prop.getProperty("ModuleLoad.Fastqc"));
				str_FastQC_Location = prop.getProperty("ModulePath.Fastqc");
				
				mTrimmomatic = Boolean.parseBoolean(prop.getProperty("ModuleLoad.Trimmomatic"));
				str_Trimmomatic_Location = prop.getProperty("ModulePath.Trimmomatic");
				
				mPicard = Boolean.parseBoolean(prop.getProperty("ModuleLoad.Picard"));
				str_Picard_Location = prop.getProperty("ModulePath.Picard");
				
				mBwa = Boolean.parseBoolean(prop.getProperty("ModuleLoad.Bwa"));
				str_BWA_Location = prop.getProperty("ModulePath.Bwa");
				
				mSamtools = Boolean.parseBoolean(prop.getProperty("ModuleLoad.Samtools"));
				str_Samtools_Location = prop.getProperty("ModulePath.Samtools");
				
				mGatk = Boolean.parseBoolean(prop.getProperty("ModuleLoad.Gatk"));
				str_Gatk_Location = prop.getProperty("ModulePath.Gatk");
				
				mSnpEff = Boolean.parseBoolean(prop.getProperty("ModuleLoad.SnpEff"));
				str_SnpEff_Location = prop.getProperty("ModulePath.SnpEff");
				
				mBedtools = Boolean.parseBoolean(prop.getProperty("ModuleLoad.Bedtools"));
				str_Bedtools_Location = prop.getProperty("ModulePath.Bedtools");
				
				str_hg38_Location= prop.getProperty("Hg38.Path");
				str_dbsnp_Location = prop.getProperty("Dbsnp.Path");
				str_hapmap_Location = prop.getProperty("Hapmap.Path");
				str_mills_Location = prop.getProperty("Mills.Path");
				str_omni_Location = prop.getProperty("Omni.Path");		
				str_phase1_Location = prop.getProperty("Phase1.Path");
				str_TruSeq3_Location = prop.getProperty("TruSeq3.Path");
				
				str_r1_filtered_trimmed_Location = prop.getProperty("ReadOne.Filtered.Trimmed.Path");
				str_r1_unfiltered_Location = prop.getProperty("ReadOne.Unpaired.Path");
				str_r2_filtered_trimmed_Location = prop.getProperty("ReadTwo.Filtered.Trimmed.Path");
				str_r2_unfiltered_Location = prop.getProperty("ReadTwo.Unpaired.Path");
			}			
			} catch (Exception e) {
				System.out.println("Exception: " + e);
			} finally {
				inputStream.close();
			}
	}
	
}
