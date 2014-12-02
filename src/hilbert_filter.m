function [ filter ] = hilbert_filter( filter_type,width )

switch filter_type
    case 'S-L'
        return
    otherwise
        if mod(width,2) ~= 0
            w = -((width-1)/2):((width-1)/2);
        else
            w = -(width/2-1):(width)/2;
        end
            filter =(1./(1i*2*pi).*sign(w))';



end

