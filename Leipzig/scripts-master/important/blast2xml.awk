{
printf("<?xml version=\"1.0\"?>\n");
printf("<BlastOutput>\n");
printf("  <BlastOutput_program>blastn</BlastOutput_program>\n");
printf("  <BlastOutput_version>BLASTN 2.2.25+</BlastOutput_version>\n");
printf("  <BlastOutput_reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>\n");
printf("  <BlastOutput_db>n/a</BlastOutput_db>\n");
printf("  <BlastOutput_query-ID>%s</BlastOutput_query-ID>\n",$1);
printf("<BlastOutput_iterations>\n");
printf("<Iteration>\n");
printf("  <Iteration_iter-num>1</Iteration_iter-num>\n");
printf("  <Iteration_query-len>%s</Iteration_query-len>\n",$8-$7);
printf("<Iteration_hits>\n");
printf("<Hit>\n");
printf("  <Hit_num>%d</Hit_num>\n",hit_num++);
printf("  <Hit_def>%s</Hit_def>\n",$2);
printf("  <Hit_len>%d</Hit_len>\n",$4);
printf("  <Hit_hsps>\n");
printf("    <Hsp>\n");
printf("      <Hsp_num>1</Hsp_num>\n");
printf("      <Hsp_bit-score>%d</Hsp_bit-score>\n",$12-$6-$7);
printf("      <Hsp_score>%d</Hsp_score>\n",$12);
printf("      <Hsp_evalue>%s</Hsp_evalue>\n",$11);
printf("      <Hsp_query-from>%s</Hsp_query-from>\n",$7);
printf("      <Hsp_query-to>%s</Hsp_query-to>\n",$8);
printf("      <Hsp_hit-from>%s</Hsp_hit-from>\n",$9);
printf("      <Hsp_hit-to>%s</Hsp_hit-to>\n",$10);
printf("      <Hsp_query-frame>???</Hsp_query-frame>\n");
printf("      <Hsp_hit-frame>??</Hsp_hit-frame>\n");
printf("      <Hsp_identity>??</Hsp_identity>\n");
printf("      <Hsp_positive>??</Hsp_positive>\n");
printf("      <Hsp_gaps>?</Hsp_gaps>\n");
printf("      <Hsp_align-len>?</Hsp_align-len>\n");
printf("    </Hsp>\n");
printf("  </Hit_hsps>\n");
printf("</Hit>\n");
printf("</Iteration_hits>\n");
printf("</Iteration>\n");
printf("</BlastOutput_iterations>\n");
printf("</BlastOutput>\n");
}
Lbil_3554       Lvir_1  98.4901 4901    58      112     20169   25069   1       4981    1e-20   4989

