for k in {0..19}
do
#    ./qfile_zaxpy WeiFreundSDML-${k}.lagrangian_multiplier.q 1.0e-3 WeiFreundSDML-${k}.diff.q WeiFreundSDML-${k}.lagrangian_multiplier.q
    ./qfile_zaxpy WeiFreundSDML-${k}.lagrangian_multiplier.q 0.0 WeiFreundSDML-${k}.diff.q --zero
done
