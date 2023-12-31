#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process make_csv {
	publishDir "${params.outdir}"
	input:
	path(fastq_input)
	output:
	path("samplelist.csv")
	
	script:
	"""
	ls -1 ${fastq_input} > sample.csv
	realpath ${fastq_input}/* > paths.csv
	paste sample.csv paths.csv > samplelist.csv
	sed -i 's/	/,/g' samplelist.csv
	sed -i '1i SampleName,SamplePath' samplelist.csv
	"""

}

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
    path("combined_ANI.csv"),emit:combined_ANI
    script:
    """
    mkdir all_assembly
    for i in ${assembly_list};do cp \$i "all_assembly"/\$i;done
    cat ${ani_list} > combined_ANI.csv
    sed -i "1i Query\ssequence\tReference\ssequence\tAverage\sNucleotide\sIdentity\tOrthologous\sMatches\tTotal\ssequence\sfragments" "combined_ANI.csv"
    """
}
process aniclustermap {
    label "medium"
    publishDir "${params.outdir}/aniclustermap",mode:"Copy"
    input:
    path (assembly_dir)
    output:
    path("ani_dir"),emit:clustermap
    path("clustermap.png"),emit:png
    script:
    """
    ANIclustermap -i ${assembly_dir} -o ani_dir
    cp ani_dir/*.png clustermap.png

    """
}
process make_report {
    label "medium"
    publishDir "${params.outdir}/reports",mode:"Copy"
    input:
    path (ANI)
    path(image)
    path(rmdfile)
    output:
    path("report.html")
    script:
    """
    cp ${rmdfile} rmdfile1.rmd
    cp ${ANI} ani.csv
    cp ${image} test.png
    Rscript -e 'rmarkdown::render(input="rmdfile1.rmd",params=list(csvfile= "ani.csv",png="test.png"),output_file="report.html")'
    """
}
workflow {
    data=Channel
	.fromPath(params.input)
    merge_fastq(make_csv(data).splitCsv(header:true).map { row-> tuple(row.SampleName,row.SamplePath)})
	if (params.trim_barcodes){
		porechop(merge_fastq.out)
		dragonflye(porechop.out,params.gsize) 
	} else {
        dragonflye(merge_fastq.out,params.gsize)           
    }
    fastANI(dragonflye.out.sample,dragonflye.out.assembly,params.reference)
    combine_outputs(dragonflye.out.assembly.collect(),fastANI.out.collect())
    aniclustermap(combine_outputs.out.all_assem)
    rmd_file=file("${baseDir}/report.rmd")
    make_report(combine_outputs.out.combined_ANI,aniclustermap.out.png,rmd_file)
}
