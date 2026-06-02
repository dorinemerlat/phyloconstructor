#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PHYLOCONSTRUCTOR } from './workflows/phyloconstructor'

workflow {
    PHYLOCONSTRUCTOR()
}