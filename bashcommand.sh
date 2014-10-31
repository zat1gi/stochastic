./astochastic | tee >(grep sign | wc -l > log1.print) >(grep sample | wc -l > log2.print) | grep : ; cat log1.print ; cat log2.print ; rm log1.print log2.print
