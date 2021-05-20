nextflow.enable.dsl=2

process HICUP {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	label 'bigMem'
	label 'multiCore'
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (hicup_args)
		val (verbose)
		val (enzyme)

	output:
	    path "*.bam",  emit: bam
		path "*.txt", emit: stats 
		path "*.html", emit: reports
		path "*.svg", emit: plots
		path "*.homer", emit: homer
		path "*.prejuicer", emit: prejuicer
		

	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:
		if (verbose){
			println ("[MODULE] HICUP ARGS: " + hicup_args)
		}

		cores = 2
		readString = ""

		// Options we add are
		hicup_options = hicup_args
		readString = reads[0] + ' ' + reads[1]

		index = params.genome["bowtie2"]
		hicup_name = name + "_" + params.genome["name"]

		//Check digest
		hicup_digest_folder = params.genome["hicup_digest"]
		myDir = file(hicup_digest_folder)
		allFiles = myDir.list()
		hicup_digest_file = ''
           
		enzyme_match_counter = 0     
		for( def file : allFiles ) {
			digest_file_enzyme = file.split('_');
			digest_file_enzyme = digest_file_enzyme[2]

			if(digest_file_enzyme == enzyme){
				enzyme_match_counter = enzyme_match_counter + 1
				hicup_digest_file = hicup_digest_folder + "/" + file
			}
		}

		if(enzyme_match_counter != 1){
			println("Matching restriction enzymes not exactly equal to 1")
			println("Specify ONE valid hi-c restriction --enzyme")
			System.exit(1)
		}

		"""
		hicup --bowtie2 bowtie2 --format sanger --index $index --threads $cores --longest 700 --shortest 50 --zip --digest $hicup_digest_file ${hicup_options} $readString
		hicup2homer *.bam
		hicup2juicer *.bam
		"""
}