package edu.rutgers.ifh;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


public class CLI {
    private String[] args = null;
    private Options options = new Options();
    private CommandLineParser parser = null;
    private CommandLine cmd = null;
    
    public String[] reads;
    public String[] ids;
    public String output;
    public String nodes;
    public String ppn;
    public String memory;
    public String wall;
    public String email;

	public CLI(String[] args) {
        this.args = args;
        /*Prints the help*/ 
        options.addOption("h", "help", false, "Show help.");
        options.addRequiredOption("o", "output", true, "Set output directory.");
        options.addRequiredOption("n", "nodes", true, "Number of nodes need it.");
        options.addRequiredOption("p", "ppn", true, "Set ppn.");   
        options.addRequiredOption("m", "memory", true, "Set required memory.");
        options.addRequiredOption("w", "wallTime", true, "Set wall time.");
        options.addRequiredOption("e", "email", true, "Set notification email.");
        
        Option R = new Option("r", "reads", true, "Set all pairs of reads.");
        R.setArgs(1000);
        R.setOptionalArg(true);
        
        Option ID = new Option("id", "sampleID", true, "Set all sample ids.");
        ID.setArgs(100);
        ID.setOptionalArg(true);
        
        options.addOption(R);
        options.addOption(ID);
    }

    public void parse() {
        parser = new DefaultParser();        
        try {
            cmd = parser.parse(options, args);

            if (cmd.hasOption("h"))
                help();

            reads = cmd.getOptionValues("r");
            ids = cmd.getOptionValues("id");
            output = cmd.getOptionValue("o");
            nodes = cmd.getOptionValue("n");
            ppn = cmd.getOptionValue("p");
            memory = cmd.getOptionValue("m");
            wall = cmd.getOptionValue("w");
            email = cmd.getOptionValue("e");
            
            if(output.charAt(output.length() -1) == '/')
            	output.substring(0, output.length() - 1);

        } catch (ParseException ex) {
        	help();
        }
    }

    private void help() {
        // This prints out some help
        HelpFormatter formater = new HelpFormatter();

        formater.printHelp("Main", options);
        System.exit(0);
    }
    
    public Options getOptions() {
		return options;
	}

	public void setOptions(Options options) {
		this.options = options;
	}

	public CommandLineParser getParser() {
		return parser;
	}

	public void setParser(CommandLineParser parser) {
		this.parser = parser;
	}

	public String[] getReads() {
		return reads;
	}

	public void setReads(String[] reads) {
		this.reads = reads;
	}

	public String[] getId() {
		return ids;
	}

	public void setId(String[] id) {
		this.ids = id;
	}

	public String getOutput() {
		return output;
	}

	public void setOutput(String output) {
		this.output = output;
	}

	public String getNodes() {
		return nodes;
	}

	public void setNodes(String nodes) {
		this.nodes = nodes;
	}

	public String getPpn() {
		return ppn;
	}

	public void setPpn(String ppn) {
		this.ppn = ppn;
	}

	public String getMemory() {
		return memory;
	}

	public void setMemory(String memory) {
		this.memory = memory;
	}

	public String getWall() {
		return wall;
	}

	public void setWall(String wall) {
		this.wall = wall;
	}

	public String getEmail() {
		return email;
	}

	public void setEmail(String email) {
		this.email = email;
	}
}