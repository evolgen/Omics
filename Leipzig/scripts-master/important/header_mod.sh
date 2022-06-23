#!/usr/bin/bash
awk '/^>/{$0=$0"_"(++i)}1' chimp_dna_sm.fa > chimp_dna_sm_mod.fa
