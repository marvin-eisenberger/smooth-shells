function blocks = blockstackmat(M,numBlocks,blockWidth)

    blocks = M * kron(ones(numBlocks,1),speye(blockWidth));

end

