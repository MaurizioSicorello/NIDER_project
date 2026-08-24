function SI_out = insert_cluster_fields(SI_out, SI_in, idx)

f = fieldnames(SI_out);

for ii = 1:numel(f)
    thisfield = f{ii};

    try
        a = SI_out.(thisfield);
        b = SI_in.(thisfield);

        if isnumeric(a) && isnumeric(b)

            if isvector(a) && numel(a) >= max(idx) && numel(b) >= 1
                SI_out.(thisfield)(idx) = b(1);

            elseif ismatrix(a) && size(a,1) >= max(idx) && size(b,1) >= 1 && size(a,2) == size(b,2)
                SI_out.(thisfield)(idx, :) = repmat(b(1, :), numel(idx), 1);
            end
        end
    catch
    end
end