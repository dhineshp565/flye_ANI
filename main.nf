#!/usr/bin/env nextflow
nextflow.enable.dsl=2


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
	porechop -i ${SamplePath} -o ${SampleName}_trimmed.fastq
	"""
}
process dragonflye {
    label "high"
    publishDir "${params.outdir}/Assembly",mode:"copy"
    input:
    tuple val(SampleName),path(SamplePath)
    val(gsize)
    output:
    val(SampleName),emit:sample
	path("${SampleName}_flye.fasta"),emit:assembly
	path("${SampleName}_flye-info.txt"),emit:flyeinfo
    script:
    """
    dragonflye --reads ${SamplePath} --outdir ${SampleName}_assembly --model r1041_e82_400bps_sup_g615 --gsize ${gsize} --nanohq --medaka 1
    mv "${SampleName}_assembly"/flye.fasta "${SampleName}"_flye.fasta
    sed -i 's/contig/${SampleName}_contig/g' "${SampleName}_flye.fasta"
    mv "${SampleName}_assembly"/flye-info.txt "${SampleName}"_flye-info.txt
    sed -i 's/contig/${SampleName}_contig/g' "${SampleName}_flye-info.txt"
    """
}

process fastANI {
    label "medium"
    publishDir "${params.outdir}/fastANI",mode:"Copy"
    input:
    val(SampleName)
    path(fasta)
    path(reference)
    output:
    path("${SampleName}_fastANI_out.txt"),emit:ANI
    script:
    """
    fastANI -q ${fasta} -r ${reference} -o ${SampleName}_fastANI_out.txt

    """
}
process combine_outputs {
    label "low"
    publishDir "${params.outdir}/combined_outputs",mode:"Copy"
    input:
    path (assembly_list)
    path(ani_list)
    output:
    path("all_assembly"),emit:all_assem
    path("combined_ANI.csv")
    script:
    """
    mkdir all_assembly
    for i in ${assembly_list};do cp \$i "all_assembly"/\$i;done
    cat ${ani_list} > combined_ANI.csv
    sed -i "1i Querysequence Reference_sequence Average_Nucleotide_Identity Orthologous_Matches Total_sequence_fragments" "combined_ANI.csv"
    """
}
process aniclustermap {
    label "medium"
    publishDir "${params.outdir}/aniclustermap",mode:"Copy"
    input:
    path (assembly_dir)
    output:
    path("ani_dir"),emit:clustermap
    script:
    """
    ANIclustermap -i ${assembly_dir} -o ani_dir
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
		dragonflye(porechop.out,params.gsize) 
	} else {
        dragonflye(merge_fastq.out,params.gsize)           
    }
    fastANI(dragonflye.out.sample,dragonflye.out.assembly,params.reference)
    combine_outputs(dragonflye.out.assembly.collect(),fastANI.out.collect())
    aniclustermap(combine_outputs.out.all_assem)
}