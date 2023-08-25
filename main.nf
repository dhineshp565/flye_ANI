#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//parse input csv file with SampleName,path to SampleName fastq

params.input='./SampleList.csv'
params.outdir='Results'
params.medaka_model=' '
params.gsize=' '
params.ANIout='sample.txt'


//merge fastq files for each SampleName and create a merged file for each SampleNames
process merge_fastq {
	publishDir "${params.outdir}/merged"
	label "low"
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path("${SampleName}.{fastq,fastq.gz}")
	
	shell:
	"""
	count=\$(ls -1 $SamplePath/*.gz 2>/dev/null | wc -l)
	
	
		if [[ "\${count}" != "0" ]];
		then
			cat $SamplePath/*.fastq.gz > ${SampleName}.fastq.gz
		
		else
			cat $SamplePath/*.fastq > ${SampleName}.fastq
		fi
	"""
}
process porechop {
	label "medium"
	publishDir "${params.outdir}/trimmed",mode:"copy",overwrite: false
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path ("${SampleName}_trimmed.fastq")
	script:
	"""
	porechop -i $SamplePath -o ${SampleName}_trimmed.fastq
	"""
}
process dragonflye {
    label "high"
    publishDir "${params.outdir}/Assembly",mode:"copy"
    input:
    tuple val(SampleName),path(SamplePath)
    val(medaka_model)
    val(gsize)
    output:
    val(SampleName),emit:SampleName
	path("${SampleName}_flye.fasta"),emit:assembly
	path("${SampleName}_flyeinfo.txt"),emit:flyeinfo
    script:
    """
    dragonflye --reads ${SampleName}.fastq --outdir ${SampleName}_assembly --model ${medaka_model} --gsize ${gsize} --nanohq --ram 4 --medaka 1
    mv "${SampleName}_assembly"/flye.fa "${SampleName}"_flye.fasta
    mv "${SampleName}_assembly"/flyeinfo.txt "${SampleName}"_flyeinfo.txt
    """
}

process fastANI {
    label "medium"
    publishDir "${params.outdir}/fastANI",mode:"Copy"
    input:
    path(query_list)
    path(reference)
    val (fastANI_out)
    output:
    path ("${fastANI_out}.txt"),emit:ANI
    script
    """
    fastANI --ql ${query_list} -r ${reference} -o ${fastANI_out}.txt
    sed '1i Query.sequence	Reference.Sequence	Average.Nucleotide.Identity	Orthologous.Matches	Total.Sequence.fragments' "${fastANI_out}.txt"
    """
}
process aniclustermap {
    label "medium"
    publishDir "${params.outdir}/aniclustermap",mode:"Copy"
    input:
    path (assembly_dir)
    path(ani_dir)
    output
    path("${fastANI_out}.png"),emit:clustermap
    script:
    """
    ANIclustermap -i ${assembly_dir} -o ${ani_dir}
    cp "${ani_dir}"/*.png ${fastANI_out}.png
    """
}
workflow {
    data=Channel
	.fromPath(params.input)
	.splitCsv(header:true)
    .map { row-> tuple(row.SampleName,row.SamplePath) }
    merge_fastq(data)

    if (params.trim_barcodes){
		porechop(merge_fastq.out)
		dragonflye(porechop.out,params.medaka_model,params.gsize)
		 
	} else {
        dragonflye(merge_fastq.out,params.medaka_model,params.gsize)
                
    }
}



