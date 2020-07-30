%Function to create a block matrices
function output = blockMatrices(matrix, rows, cols)
    a  = size(matrix, 1);
    b  = size(matrix, 2);
    c = floor(a/rows);
    d = rem(a, rows);
    partition_a = ones(1, rows)*c;
    partition_a(1:d) = partition_a(1:d)+1;
    e = floor(b/cols);
    f = rem(b, cols);
    partition_b = ones(1, cols)*e;
    partition_b(1:f) = partition_b(1:f)+1;
    %output = ones(rows, cols, partition_a, partition_b);
    output = mat2cell(matrix, partition_a, partition_b);   
end