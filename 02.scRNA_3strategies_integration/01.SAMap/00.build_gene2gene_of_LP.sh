
file1='00.data0/lamprey.fa'
type1='nucl'
id1='LP'

######################################################################################
file2='00.data0/freshwater_butterflyfish.fa'
type2='nucl'
id2='BF'
bash ./00.basic_BLASTP_one2one.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2


file2='00.data0/indian_medaka.fa'
type2='nucl'
id2='Om'
bash ./00.basic_BLASTP_one2one.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2


file2='00.data0/tarpon.fa'
type2='nucl'
id2='Ma'
bash ./00.basic_BLASTP_one2one.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2

file2='00.data0/whitespotted_bambooshark.fa'
type2='nucl'
id2='Wb'

bash ./00.basic_BLASTP_one2one.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2

file2='00.data0/lungfish.fa'
type2='nucl'
id2='LF'

bash ./00.basic_BLASTP_one2one.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2

file2='00.data0/mississippi_paddlefish.fa'
type2='nucl'
id2='Pd'

bash ./00.basic_BLASTP_one2one.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2

file2='00.data0/gray_bichir.fa'
type2='nucl'
id2='Br'

bash ./00.basic_BLASTP_one2one.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2

######################################################################################
file2='00.data0/human.fa'
type2='nucl'
id2='Hm'
bash ./00.basic_BLASTP_one2one.sh --tr1 $file1 --t1 $type1 --n1 $id1 --tr2 $file2 --t2 $type2 --n2 $id2
