/* Exercise publication with the same Nextflow 23.10.1 runtime as Crest. */
nextflow.enable.dsl = 2

process PUBLICATION_CHECK {
    publishDir "${launchDir}/publication-check", mode: 'symlink', overwrite: true
    output:
    path('publication-file.txt')
    path('publication-directory')
    script:
    '''
    printf 'publication check\n' > publication-file.txt
    mkdir publication-directory
    for i in $(seq 1 1000); do
        printf '%s\n' "$i" > "publication-directory/$i.txt"
    done
    '''
}

workflow {
    PUBLICATION_CHECK()
}

workflow.onComplete {
    if (workflow.success) {
        def root = file("${launchDir}/publication-check")
        assert java.nio.file.Files.isSymbolicLink(root.resolve('publication-file.txt'))
        assert java.nio.file.Files.isSymbolicLink(root.resolve('publication-directory'))
        assert root.resolve('publication-directory/1000.txt').text.trim() == '1000'
        println 'PASS: file and 1,000-file directory published as symbolic links'
    }
}
