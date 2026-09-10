# RFQSC: copia local del elemento RFQ de TRANSOPTR

Fuente: `/home/cvalerio/work1/transoptr/transoptr-master/src/RFQ.f`.
Autor identificado en la fuente: Olivier Shelbaya, TRIUMF.
Se conserva la licencia GPL de la distribución en `LICENSE.txt`.

Modificaciones: 9 de septiembre de 2026. Rutinas RFQ, RFQn y SCRFQ renombradas
RFQSC, RFQSCn y SCRFQSC para convivir con los objetos instalados. En SCRFQSC se
sustituye BNCHARGE, una variable local sin asignación, por CURRENT de COMMON/MOM/.
En modo 5 CURRENT contiene carga de bunch en coulombs. La fórmula QSC resultante
es equivalente a la de SC.f en modo 5. El hash de la fuente original está en la
cabecera y el comparador registra el hash de esta copia en comparison.json.

La copia se agrega al sy.f generado solo cuando `transoptr.rfq_source_override`
está configurado. El comparador de corriente lo activa en ambos casos. No se
modificaron archivos ni objetos de la instalación global de TRANSOPTR.
