package edu.rutgers.ifh;

import org.apache.commons.io.FilenameUtils;

public class QC_WES {

	private QC_Tools _tools;
	private CLI _cli;
	private String R1 = "";
    private String _dbUsername = "";
    private String _dbPassword = "";

	QC_WES(QC_Tools tools, CLI cli, String user, String password) {
		_tools = tools;
		_cli = cli;
		_dbUsername = user;
		_dbPassword = password;
	}

	String CreateHeader() {

		String str_Bin_Bash = "#!/bin/bash \n \n";

		String str_Partition = "#SBATCH --partition=main\n";
		String str_Reque = "#SBATCH --requeue \n";
		String str_JobName = "#SBATCH --job-name=" + _tools.str_job_name + "\n";
		String str_Nodes = "#SBATCH --nodes=" + _cli.nodes + "\n";
		String str_Tasks = "#SBATCH --ntasks=" + "1 \n";
		String str_Cpu = "#SBATCH --cpus-per-task=" + "1\n";
		String str_Memory = "#SBATCH --mem=" + _cli.memory + "\n";
		String str_Time = "#SBATCH --time=" + _cli.wall + "\n";
		String str_Output = "#SBATCH --output=slurm.%N.%j.out \n";
		String str_Error = "#SBATCH --error=slurm.%N.%j.err \n";
		String str_Email = "#SBATCH --mail-user=" + _cli.email + "\n";
		String str_Export = "#SBATCH --export=ALL \n";

		String str_Header = str_Bin_Bash + str_Partition + str_Reque + str_JobName + str_Nodes + str_Tasks + str_Cpu
				+ str_Memory + str_Time + str_Output + str_Error + str_Email + str_Export;

		String str_ModulePurge = "module purge";
		
		String str_ModuleCommunity = "module use /projects/community/modulefiles";
		
		String str_Modules = "module load ";

		if (_tools.mFastqc)
			str_Modules += _tools.str_FastQC_Location + " ";
		if (_tools.mTrimmomatic)
			str_Modules += _tools.str_Trimmomatic_Location + " ";
		if (_tools.mPicard)
			str_Modules += _tools.str_Picard_Location + " ";
		if (_tools.mBwa)
			str_Modules += _tools.str_BWA_Location + " ";
		if (_tools.mSamtools)
			str_Modules += _tools.str_Samtools_Location + " ";
		if (_tools.mGatk)
			str_Modules += _tools.str_Gatk_Location + " ";
		if (_tools.mSnpEff)
			str_Modules += _tools.str_SnpEff_Location + " ";
		if (_tools.mBedtools)
			str_Modules += _tools.str_Bedtools_Location + " ";
				
		return str_Header + "\n" + str_ModulePurge + "\n" + str_ModuleCommunity + "\n" + str_Modules +"\n"; 
	}

	String CreateDir(int id) {
		String command = "mkdir -p " + _cli.output + "/" + _cli.ids[id] + "/ \n";
		command += "mkdir -p " + _cli.output + "/" + _cli.ids[id] + "/indels \n";
		command += "mkdir -p " + _cli.output + "/" + _cli.ids[id] + "/snps \n\n";
		return command;
	}

	String CdDir(int id) {
		String command = "cd " + _cli.ids[id] + "/ \n";
		return command;
	}

	String QualityControl(int r1, int r2, int id) {
		String command = "srun ";

		if (_tools.mFastqc)
			command += "fastqc ";
		else
			command += _tools.str_FastQC_Location;

		return command + " --quiet -o " + ". " + _cli.reads[r1] + " "
				+ _cli.reads[r2] + "\n";
	}

	String Trimmomatic(int r1, int r2, int id) {
		String command = "srun ";

		if (_tools.mTrimmomatic)
			command += "trimmomatic";
		else
			command += "java -jar " + _tools.str_Trimmomatic_Location;

		return command + " PE " + _cli.reads[r1] + " " + _cli.reads[r2] + " " + R1 + "_1_filtered_trimmed.fq.gz"
				+ " " + R1 + "_1_unpaired.fq.gz" + " " + R1 + "_2_filtered_trimmed.fq.gz" + " "
				+ R1 + "_2_unpaired.fq.gz"
				+ " ILLUMINACLIP:" + _tools.str_TruSeq3_Location + ":2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:6\n";
	}

	String Indexing() {
		String command = "srun ";

		if (_tools.mBwa)
			command += "bwa ";
		else
			command += _tools.str_BWA_Location;

		return command + " index " + _tools.str_hg38_Location + "\n";
	}

	String ALignmentOne(int id) {
		String command = "srun ";

		if (_tools.mBwa)
			command += "bwa";
		else
			command += _tools.str_Bedtools_Location;

		//+ _cli.output + "/" + _cli.ids[id] + "/"
		return command + " mem -M -t 10 " + _tools.str_hg38_Location + " " + R1 + "_1_filtered_trimmed.fq.gz"
				+ " " +  R1 + "_2_filtered_trimmed.fq.gz" + " > "  + R1
				+ "_bwa.sam" + "\n";
	}

	String ALignmentTwo(int id) {
		String command = "srun ";

		if (_tools.mBwa)
			command += "bwa";
		else
			command += _tools.str_BWA_Location;

		return command + " mem -R '@RG\\tID:dip1\\tSM:dip1\\tPL:ILLUMINA' " + _tools.str_hg38_Location
				+ R1 + "_1_filtered_trimmed.fq.gz" + " " + R1 + "_2_filtered_trimmed.fq.gz" + " > " + R1 + "_bwa.sam" + "\n";
	}

	String SortBam(int id) {

		String command = "srun ";

		if (_tools.mPicard)
			command += "picard";
		else
			command += "java -Xmx10g -jar " + _tools.str_Picard_Location;

		return command + " SortSam VALIDATION_STRINGENCY=SILENT I=" + R1
				+ "_bwa.sam" + " O=" +  R1 + "_Sort_bwa.sam SORT_ORDER=coordinate" + "\n";
	}

	String RemoveDuplicates(int id) {
		String command = "srun ";
		if (_tools.mPicard)
			command += "picard ";
		else
			command += "java -Xmx10g -jar " + _tools.str_Picard_Location;

		return command + " MarkDuplicates VALIDATION_STRINGENCY=SILENT I=" + R1
				+ "_Sort_bwa.sam " + " O=" + R1
				+ "_PCR_bwa.sam REMOVE_DUPLICATES=true M=" + R1
				+ "_pcr_bwa.metrics" + "\n";
	}

	String CollectInsertSizeMetrics(int id) {
		String command = "srun ";
		if (_tools.mPicard)
			command += "picard ";
		else
			command += "java -Xms1g -Xmx4g -jar " + _tools.str_Picard_Location;

		return command + " CollectInsertSizeMetrics METRIC_ACCUMULATION_LEVEL=ALL_READS OUTPUT=" + R1 + "_rmdup_insertSize.txt HISTOGRAM_FILE=" +
				R1 + "_rmdup_metrics.pdf INPUT=" + R1 + "_PCR_bwa.sam"
				+ "\n";
	}

	String SamtoolsSort(int id) {
		String command = "srun ";
		if (_tools.mSamtools)
			command += "samtools ";
		else
			command += "java -Xms1g -Xmx4g -jar " + _tools.str_Samtools_Location;

		return command + " sort " + R1 + "_PCR_bwa.sam -o " + R1 + "_PCR_bwa_sorted.bam" + "\n";
	}

	String SamtoolsIndex(int id) {
		String command = "srun ";
		if (_tools.mSamtools)
			command += "samtools ";
		else
			command += "java -Xms1g -Xmx4g -jar " + _tools.str_Samtools_Location;

		return command + " index " + R1 + "_PCR_bwa_sorted.bam" + "\n";
	}

	String CreateRealignmentTargets(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T RealignerTargetCreator -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T RealignerTargetCreator -R ";

		return command + _tools.str_hg38_Location + " -I " + R1
				+ "_PCR_bwa_sorted.bam -o realignment_targets.list " + "\n";
	}

	String RealignIndels(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T IndelRealigner -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T IndelRealigner -R ";

		return command + _tools.str_hg38_Location + " -I " + R1
				+ "_PCR_bwa_sorted.bam -targetIntervals realignment_targets.list -o realigned_reads.bam" + "\n";

	}

	String BaseQualityScoreRecalibrationBQSR_1(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T BaseRecalibrator -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T BaseRecalibrator -R ";

		return command + _tools.str_hg38_Location + " -I realigned_reads.bam -knownSites " + _tools.str_dbsnp_Location + " -knownSites "
				+ _tools.str_hapmap_Location + " -knownSites " + _tools.str_mills_Location + " -o recal_data.table" + "\n";
	}

	String BaseQualityScoreRecalibrationBQSR_2(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T BaseRecalibrator -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T BaseRecalibrator -R ";

		return command + _tools.str_hg38_Location + " -I realigned_reads.bam -knownSites " + _tools.str_dbsnp_Location + " -knownSites "
				+ _tools.str_hapmap_Location + " -knownSites " + _tools.str_mills_Location + " -BQSR recal_data.table -o post_recal_data.table" + "\n";
	}

	String AnalyzeCovariates(int id) {

		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T AnalyzeCovariates -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T AnalyzeCovariates -R ";

		return command + _tools.str_hg38_Location + " -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf\n";
	}

	String ApplyBQSR(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T PrintReads -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T PrintReads -R ";

		return command + _tools.str_hg38_Location + " -I realigned_reads.bam -BQSR recal_data.table -o recal_reads.bam" + "\n";
	}

	String CallVariants(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T HaplotypeCaller -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T HaplotypeCaller -R ";

		return command + _tools.str_hg38_Location + " -I recal_reads.bam -o raw_variants_recal.vcf" + "\n";
	}

	String ExtractRawSNPs(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T SelectVariants -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T SelectVariants -R ";

		return command + _tools.str_hg38_Location + " --variant raw_variants_recal.vcf -selectType SNP -o raw_snps_recal.vcf" + "\n";
	}

	String ExtractRawIndels(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T SelectVariants -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T SelectVariants -R ";

		return command + _tools.str_hg38_Location + " --variant raw_variants_recal.vcf -selectType INDEL -o raw_indels_recal.vcf" + "\n";
	}

	String SNPRecalibrationModel(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T VariantRecalibrator -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T VariantRecalibrator -R ";

		return command + _tools.str_hg38_Location + "--input raw_snps_recal.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 "
				+ _tools.str_hapmap_Location + " -resource:omni,known=false,training=true,truth=true,prior=12.0 "
				+ _tools.str_omni_Location + " -resource:1000G,known=false,training=true,truth=false,prior=10.0 "
				+ _tools.str_phase1_Location + " -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "
				+ _tools.str_dbsnp_Location
				+ " -an QD -an FS -an MQRankSum -an ReadPosRankSum -an SOR -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile "
				+ "recalibrate_SNP.recal -tranchesFile "
				+ "recalibrate_SNP.tranches --maxGaussians 4 -rscriptFile "
				+ "recalibrate_SNP_plots.R" + "\n";
	}

	String ApplyRecalibrationSNPs(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T ApplyRecalibration -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T ApplyRecalibration -R ";

		return command + _tools.str_hg38_Location + "--input "
				+ "raw_snps_recal.vcf -mode SNP -ts_filter_level 99.0 -recalFile " 
				+ "recalibrate_SNP.recal -tranchesFile " 
				+ "recalibrate_SNP.tranches -o " 
				+ "recalibrated_snps_raw_indels.vcf" + "\n";
	}

	String IndelRecalibrationModel(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T VariantRecalibrator -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T VariantRecalibrator -R ";

		return command + _tools.str_hg38_Location + " --input "
				+ "raw_indels_recal.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "
				+ _tools.str_dbsnp_Location + " -resource:mills,known=false,training=true,truth=true,prior=12.0 "
				+ _tools.str_mills_Location
				+ " -an QD -an FS -an SOR -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0  -recalFile "
				+ "recalibrate_INDEL.recal -tranchesFile "
				+ "recalibrate_INDEL.tranches --maxGaussians 4 -rscriptFile "
				+ "recalibrate_INDEL_plots.R" + "\n";
	}

	String ApplyRecalibrationIndels(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T ApplyRecalibration -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T ApplyRecalibration -R ";
		return command + _tools.str_hg38_Location + " --input "
				+ "raw_indels_recal.vcf -mode INDEL -ts_filter_level 99.0 -recalFile "
				+ "recalibrate_INDEL.recal -tranchesFile "
				+ "recalibrate_INDEL.tranches -o recalibrated_variants.vcf"
				+ "\n";
	}

	String ExtractFilteredSNPs(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T SelectVariants -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T SelectVariants -R ";

		return command + _tools.str_hg38_Location + " --variant "
				+ "recalibrated_snps_raw_indels.vcf -selectType SNP -o "
				+ "filtered_snps_final.vcf" + "\n";
	}

	String ExtractFilteredIndels(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T SelectVariants -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T SelectVariants -R ";

		return command + _tools.str_hg38_Location + " --variant "
				+ "recalibrated_variants.vcf -selectType INDEL -o "
				+ "filtered_indels_final.vcf" + "\n";
	}

	String EvaluateFilteredSNVsAgainstSNVDatabase(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T VariantEval -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T VariantEval -R ";

		return command + _tools.str_hg38_Location + " -eval "
				+ "recalibrated_variants.vcf -D " + _tools.str_dbsnp_Location
				+ " -noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants -EV MultiallelicSummary -o "
				+ "Eval_report.txt" + "\n";
	}

	String TabulateAnnotationValuesSNPs(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T VariantsToTable -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T VariantsToTable -R ";

		return command + _tools.str_hg38_Location + " --variant "
				+ "filtered_snps_final.vcf -F CHROM -F POS -F ID -F ALT -F QUAL -F BaseQRankSum -F DP -F DP_Orig -F FS -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR -F set -F varType -GF GQ -GF DP -o "
				+ "combinedSNPeval.txt" + "\n";
	}

	String TabulateAnnotationValuesIndels(int id) {
		String command = "srun ";

		if (_tools.mGatk)
			command += "gatk -T VariantsToTable -R ";
		else
			command += "java -jar " + _tools.str_Gatk_Location + " -T VariantsToTable -R ";

		return command + _tools.str_hg38_Location + " --variant "
				+ "filtered_indels_final.vcf -F CHROM -F POS -F ID -F ALT -F QUAL -F BaseQRankSum -F DP -F DP_Orig -F FS -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR -F set -F varType -GF GQ -GF DP -o "
				+ "combinedIndeleval.txt" + "\n";
	}

	String ComputeCoverageStatistics(int id) {
		String command = "srun ";

		if (_tools.mBedtools)
			command += "bedtools ";
		else
			command += _tools.str_Bedtools_Location;

		return command + " genomecov -bga -ibam recal_reads.bam > "
				+ "genomecov.bedgraph" + "\n";
	}

	String AnnotateSNPsPredictEffects(int id) {
		String command = "srun ";

		if (_tools.mSnpEff)
			command += "snpEff ";
		else
			command += "java -jar " + _tools.str_SnpEff_Location;

		return command + " -v GRCh38.86 filtered_snps_final.vcf > "
				+ "snps/filtered_snps_final.ann.vcf" + "\n";
	}

	String AnnotateIndelsPredictEffects(int id) {
		String command = "srun ";

		if (_tools.mSnpEff)
			command += "snpEff ";
		else
			command += "java -jar " + _tools.str_SnpEff_Location;

		return command + " -v GRCh38.86 filtered_indels_final.vcf > "
				+ "indels/filtered_indels_final.ann.vcf" + "\n";

	}
	
	String Vcf2Sql(String filename, int id ) {
		String command = "srun java -jar " + _tools.str_Vcf2Sql_Location + " " +_tools.str_DB_Host_Name + " " + _tools.int_DB_port + " " + _tools.str_DB_Table + " " + _dbUsername + " " + _dbPassword + " " + filename + " " + id + " snpEff false\n";
		return command;
	}

	String RunWesPipeline(int r1, int r2, int id) {
		String fullCommand = "";

		R1 = FilenameUtils.getBaseName(_cli.reads[r1]);
		R1 = R1.split("\\.")[0];
		
		fullCommand += CreateHeader();
		
		fullCommand += QualityControl(r1,r2,id);
		 
		fullCommand += Trimmomatic(r1, r2, id);
		 
		fullCommand += Indexing();
		  
		fullCommand += ALignmentOne(id);
		  
		fullCommand += ALignmentTwo(id);
		 
		fullCommand += SortBam(id);
		  
		fullCommand += RemoveDuplicates(id);
		 
		fullCommand += CollectInsertSizeMetrics(id);
		 
		fullCommand += SamtoolsSort(id);
		 
		fullCommand += SamtoolsIndex(id);
		  
		fullCommand += CreateRealignmentTargets(id);
		  
		fullCommand += RealignIndels(id);
		  
		fullCommand += BaseQualityScoreRecalibrationBQSR_1(id);
		  
	    fullCommand += BaseQualityScoreRecalibrationBQSR_2(id);
		  
		fullCommand += AnalyzeCovariates(id);
		  
		fullCommand += ApplyBQSR(id);
		  
		fullCommand += CallVariants(id);
		  
		fullCommand += ExtractRawSNPs(id);
		  
		fullCommand += ExtractRawIndels(id);
		  
		fullCommand += ExtractRawIndels(id);
		  
		fullCommand += SNPRecalibrationModel(id);
		  
		fullCommand += ApplyRecalibrationSNPs(id);
		
		fullCommand += IndelRecalibrationModel(id);
		  
		fullCommand += ApplyRecalibrationIndels(id);
		  
		fullCommand += ExtractFilteredSNPs(id);
		 
		fullCommand += ExtractFilteredIndels(id);
		 
		fullCommand += EvaluateFilteredSNVsAgainstSNVDatabase(id);
		 
		fullCommand += TabulateAnnotationValuesSNPs(id);
		  
		fullCommand += TabulateAnnotationValuesIndels(id);
		  
		fullCommand += ComputeCoverageStatistics(id);
		  
		fullCommand += AnnotateSNPsPredictEffects(id);
		 
		fullCommand += AnnotateIndelsPredictEffects(id);
		
		fullCommand += Vcf2Sql("snps/filtered_snps_final.ann.vcf", id);
		
		fullCommand += Vcf2Sql("indels/filtered_indels_final.ann.vcf", id);
		 
		fullCommand += "\nsleep 10\n" + "sacct --format=MaxRSS,MaxDiskRead,MaxDiskWrite,Elapsed,Nodelist -j $SLURM_JOBID\n" + "sleep 2\n";
		  
		return fullCommand;
	}
}
