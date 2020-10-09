package edu.rutgers.ifh;

import java.io.ByteArrayInputStream;
import java.io.Console;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import com.jcraft.jsch.Channel;
import com.jcraft.jsch.ChannelSftp;
import com.jcraft.jsch.JSch;
import com.jcraft.jsch.JSchException;
import com.jcraft.jsch.Session;
import com.jcraft.jsch.SftpException;

public class Main {
	
	public static void main (String[] args) {
		
		PrintStream commander = null;
		Channel obj_Channel_sftp = null;
		Channel obj_Channel_ssh = null;
		OutputStream inputstream_for_the_channel = null;
		Session obj_JSch_Session = null;
        CLI cli = new CLI(args);
        QC_Tools tools = null;
		JSch obj_JSch = null;
		ChannelSftp sftp = null;
		
        String hpcUsername = "";
        String hpcPassword = "";
        String dbUsername = "";
        String dbPassword = "";
        
		try {
			cli.parse();			
			tools = new QC_Tools();
			tools.LoadPropertiesFiles();
			obj_JSch = new JSch(); 
			
	        Console console = System.console();
	        
	        hpcUsername = console.readLine("HPC Username: ");
	        hpcPassword = new String(console.readPassword("HPC Password: "));
	        dbUsername = console.readLine("DB Username: ");
	        dbPassword = new String(console.readPassword("DB Password: "));
	        
	        	
	        obj_JSch_Session = obj_JSch.getSession(hpcUsername, tools.str_Cluster_Host_Name, tools.int_Cluster_Host_port);
	        
	        obj_JSch_Session.setPassword(hpcPassword);
	        obj_JSch_Session.setConfig("StrictHostKeyChecking", "no");	        
	        obj_JSch_Session.connect(10 * 10000);	// 10 Seconds time out
        
	        //Create Channel
	        obj_Channel_sftp = obj_JSch_Session.openChannel("sftp");
	        obj_Channel_ssh = obj_JSch_Session.openChannel("shell");
	        
	        
	        inputstream_for_the_channel = obj_Channel_ssh.getOutputStream();
	        commander = new PrintStream(inputstream_for_the_channel, true);
	        obj_Channel_ssh.setOutputStream(System.out, true);
	        
	        obj_Channel_ssh.connect();
	        obj_Channel_sftp.connect(); 
	        
	        sftp = (ChannelSftp) obj_Channel_sftp;
	       
        QC_WES qWes = new QC_WES(tools, cli, dbUsername, dbPassword);
        
        int nReads = cli.reads.length / 2;
        int nIds = cli.ids.length;
                
        if(nReads != nIds) {
        	System.out.println("Missmath number of reads and ids");
        	System.exit(0);
        }
        
        int R1 = 0;
        int R2 = 1;
        int ID = 0;        
		
		try {
			sftp.mkdir(cli.output);
		} catch (SftpException e) {
			System.out.println("Directory already exsists..");
		} 
		
		sftp.cd(cli.output);
		  
		for(int i =0; i < nReads; i++) { 
			
			try {
				sftp.mkdir(cli.ids[ID]);
				sftp.mkdir(cli.ids[ID]+"/snps");
				sftp.mkdir(cli.ids[ID]+"/indels");
			} catch (SftpException e) {
				// TODO Auto-generated catch block
				System.out.println("Directory already exsists..");
			}
			
			sftp.cd(cli.ids[ID]); 
			
			try {
				sftp.rm(tools.str_job_name + ".sh");
			} catch (SftpException e) {
				System.out.println("The script does not exists..");
			}
				
			sftp.put (new ByteArrayInputStream (qWes.RunWesPipeline(R1, R2, Integer.parseInt(cli.ids[ID])).getBytes()), tools.str_job_name + ".sh" );
			sftp.chmod(775, tools.str_job_name + ".sh"); 
			sftp.cd(".."); 			
			R1++; R2++; ID++;
		}
		
		Thread.sleep(10000);
		
		commander.println("cd " + cli.output + "/");
		 
		for(int i =0; i < nReads; i++) { 
			commander.println(qWes.CdDir(i));
			commander.println("sbatch " + tools.str_job_name + ".sh");
			commander.println("cd .."); 
			commander.println("echo $PATH"); 
		}
        
        Thread.sleep(10000);
        
        obj_Channel_sftp.disconnect();
        obj_Channel_ssh.disconnect();
        obj_JSch_Session.disconnect();
        
		} catch (JSchException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SftpException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
