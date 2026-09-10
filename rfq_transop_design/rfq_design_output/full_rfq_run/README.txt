Full RFQ TRANSOPTR run folder
=============================

Files in this folder:
  data.dat  - TRANSOPTR input deck/template for the complete RFQ run
  sy.f      - generated RFQ call with the final point count and length
  fort.75   - generated RFQ coefficient table for the complete RFQ

Cells: 110
Final z: 2.7841462609E+02 cm
Final energy: 1.1941275801E+00 MeV
data.dat source: not configured
optr source: not copied

Typical use:
  cd /home/cvalerio/work1/machinlearning/rfq_full_sequence/rfq_transop_design/rfq_design_output/full_rfq_run
  ./optr

If data.dat only contains WARNING comments, set transoptr.template_data
in the design YAML and rerun the designer.
