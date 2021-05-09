nextflow.enable.dsl=2

process HICUP {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	label 'bigMem'
	label 'multiCore'
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		//val (bowtie2_args)
		val (hicup_args)
		val (verbose)
		val (enzyme)

	output:
	    path "*.bam",  emit: bam
		path "*.txt", emit: stats 
		//path "*stats.txt", emit: stats 
		path "*.html", emit: reports
		path "*.svg", emit: plots



	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:
		if (verbose){
			//println ("[MODULE] BOWTIE2 ARGS: " + bowtie2_args)
			println ("[MODULE] HICUP ARGS: " + hicup_args)
		}

		//cores = 8
		cores = 2
		readString = ""

		// Options we add are
		// bowtie_options = bowtie2_args
		hicup_options = hicup_args
		// bowtie_options +=  " --no-unal " // We don't need unaligned reads in the BAM file
		//hicup_options += "-c hicup.conf"
		
		//if (reads instanceof List) {
		//	readString = "-1 " + reads[0] + " -2 " + reads[1]
		//	bowtie_options += " --no-discordant --no-mixed " // just output properly paired reads
		//}
		//else {
		//	readString = "-U " + reads
		//}

		readString = reads[0] + ' ' + reads[1]
		println('READSTRING ' + readString)

		index = params.genome["bowtie2"]
		//bowtie_name = name + "_" + params.genome["name"]
		hicup_name = name + "_" + params.genome["name"]

		//Check digest
		hicup_digest_folder = params.genome["hicup_digest"]
		myDir = file(hicup_digest_folder)
		println hicup_digest_folder
		println enzyme
		allFiles = myDir.list()
		hicup_digest_file = ''
           
		enzyme_match_counter = 0     
		for( def file : allFiles ) {
			println file
			digest_file_enzyme = file.split('_');
			digest_file_enzyme = digest_file_enzyme[2]

			if(digest_file_enzyme == enzyme){
				enzyme_match_counter = enzyme_match_counter + 1
				hicup_digest_file = hicup_digest_folder + "/" + file
			}
		}

		if(enzyme_match_counter != 1){
			println("Matching restriction enzymes not exactly equal to 1")
			System.exit(1)
		}




		"""
		#  mls
		
		#module load bowtie2
		#  module load samtools
		
		hicup --bowtie2 bowtie2 --format sanger --index $index --threads $cores --longest 700 --shortest 50 --zip --digest $hicup_digest_file ${hicup_options} $readString

		"""

}