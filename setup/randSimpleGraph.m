function A = randSimpleGraph(n)
%RANDSIMPLEGRAPH Generate a random simple graph.
    while 1
        A = round(rand(n));
        A = triu(A) + triu(A,1)';
        A = A - diag(diag(A));
        if (calc_lambda(A) > 0.0001)
            break;
        end
    end
end