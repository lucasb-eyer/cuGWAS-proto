test=test1 
dir=${HOME}/rt203005_FGLS/trunk/
tc_dir=${HOME}/rt203005_FGLS/trunk/${test}
n=2000
p=4
h=.2
m=1000
t=20

echo "${m} ${t} 250 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e m ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 500 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e m ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 1000 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e m ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 250 20 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e m ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 500 20 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e m ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 1000 20 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e m ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 250 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e t ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 500 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e t ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 1000 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e t ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 250 20 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e t ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 500 20 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e t ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 1000 20 ${n} ${p} ${h}" |  ${dir}/bin/driver.x e t ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 250 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x c m ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 500 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x c m ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 1000 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x c m ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 250 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x c t ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 500 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x c t ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out
echo "${m} ${t} 1000 1 ${n} ${p} ${h}" |  ${dir}/bin/driver.x c t ${tc_dir}/X.in ${tc_dir}/Y.in ${tc_dir}/Phi.in ${tc_dir}/b.out >> ${HOME}/out/test_1.out



