function [U,SIGMA,V,time] = PictureSVD(picture,iter,full,rank)
    picture = imread(picture);
    picture = rgb2gray(picture);
    picture = im2double(picture);
    %imwrite(picture,'NewOriginal.jpg')
    %title('Original Picture');
    tic
    [U,SIGMA,V] = TruncatedSVD(picture,iter,full,rank);
    time = toc;
    picture = U*SIGMA*V';
    name = [num2str(iter) 'Iterations' num2str(full) 'Full' num2str(rank) 'Rank' num2str(time) 'sec(Untruncated).jpg'];
    imwrite(picture,name)
    %title(['Iterations =',num2str(iter),', Full =',num2str(full),', Rank =',num2str(rank),', Time =',num2str(time),'(Untruncated).jpg'])
    S = SIGMA;
    S((rank+1):end,:)=0;
    S(:,(rank+1):end)=0;
    picture = U*S*V';
    name = [num2str(iter) 'Iterations' num2str(full) 'Full' num2str(rank) 'Rank' num2str(time) 'sec(Truncated).jpg'];
    imwrite(picture,name)
    %title(['Iterations =',num2str(iter),', Full =',num2str(full),', Rank =',num2str(rank),', Time =',num2str(time),'(Truncated).jpg'])
end