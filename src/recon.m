function recon = recon(proj)
%RECON 2D reconstruction
%   
[width,projnum] = size(proj);
filter_type = 'R-L';

[ filter ] = hilbert_filter( filter_type,width );

p2_f = zeros(width,projnum);
for i = 1:projnum
    p2_f(:,i) =fft(proj(:,i)) .* filter;
end

agnle_range = linspace(0,projnum*180/(projnum+1),projnum);

p2_f = real(ifft(p2_f));
recon=iradon(p2_f,agnle_range,'none');%

end

