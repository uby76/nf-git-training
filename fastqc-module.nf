process FASTQC {
    publishDir "results/fastqc", mode: "copy"

    input:
        tuple val(ID), path(reads)

    output:
        path("${ID}_fastqc")

    script:
    """
    mkdir ${ID}_fastqc
    fastqc ${reads} -o ${ID}_fastqc
    """
}
