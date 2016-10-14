%CT - Any sp3 hybridised Carbon
%HC - Any aliphatic Hydrogen
%CT-HC=0.10900
leng=14;
bread=16;
lconst=0.4943;
natoms=11;
nrows1 = leng*bread;
nrows2 = (natoms+1)*nrows1;
nrows = (natoms+4)*nrows1;
xx=zeros(nrows,3);
offsetz=5
p=1;

%Odd layers
for i=0:sqrt(3):(bread/2-1)*sqrt(3)
    for j=1:leng
        xx(p,:)=[j*lconst i*lconst offsetz];  %Sulphur
        cr1=0;
        for (k=1:2:(natoms))
            xx(p+k,:)=[(j-0.08387)*lconst (i+cr1*0.0646)*lconst...
                offsetz+k*0.12916];
            cr1=cr1+1;
        end
        cr2=0;
        for (l=2:2:natoms)
            xx(p+l,:)=[j*lconst (i+cr2*0.0646)*lconst... 
                offsetz+l*0.12916];
            cr2=cr2+1;
        end
        p=p+natoms+1;
    end
    p=(p+(natoms+1)*(leng));
end
p=(natoms+1)*leng+1;

%Even Layers
for i=(sqrt(3)/2):sqrt(3):(bread-1)*sqrt(3)/2
    %i=in*(sqrt(3)/2);
    for j=1:leng
        xx(p,:)=[(j+0.5)*lconst i*lconst offsetz]; %Sulphur
        cr3=0;
        for (m=1:2:natoms)
            xx(p+m,:)=[(j+0.5-0.08387)*lconst (i+cr3*0.0646)*lconst... 
                offsetz+m*0.12916];
            cr3=cr3+1;
        end
        cr4=0;
        for (n=2:2:natoms)
            xx(p+n,:)=[(j+0.5)*lconst (i+cr4*0.0646)*lconst... 
                offsetz+n*0.12916];
            cr4=cr4+1
        end
        p=p+natoms+1;
    end
    p=p+(natoms+1)*leng;
end

ctr2=1;
%Putting in the Hydrogen atoms!!
for t=(natoms+1):(natoms+1):nrows2
    xx(nrows2+ctr2,:) = [xx(t,1) xx(t,2)+0.077 xx(t,3)+0.077];
    xx(nrows2+ctr2+1,:) = [xx(t,1)-0.067 xx(t,2)-0.0386 xx(t,3)+0.077];
    xx(nrows2+ctr2+2,:) = [xx(t,1)+0.067 xx(t,2)-0.0386 xx(t,3)+0.077];
    ctr2=ctr2+3; 
end


apratim = strcat('samr','.gro');
fid = fopen(apratim,'W');
fprintf(fid,'%s\n', 'SAM using Methyl head');
fprintf(fid,'%d\n',nrows);
ctr=0;
ctr3=0;
ctr4=0;

for ap=1:(nrows2/(natoms+1))
    ctr=ctr + 1;
    ctr4=ctr4+1;
    fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
      ,'RS','S',ctr4,xx(ctr,1),xx(ctr,2),xx(ctr,3));
  ctr5=0;
    for c=ap+1:ap+natoms
        ctr=ctr+1;
        ctr4=ctr4+1;
        ctr5=ctr5+1;
        fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
        ,'RS',strcat('C',num2str(ctr5)),ctr4,xx(ctr,1),xx(ctr,2),xx(ctr,3));
    end
    
  fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
      ,'RS','H',ctr4+1,xx(nrows2+ctr3+1,1),...
      xx(nrows2+ctr3+1,2),xx(nrows2+ctr3+1,3));
  fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
      ,'RS','H',ctr4+2,xx(nrows2+ctr3+2,1),...
      xx(nrows2+ctr3+2,2),xx(nrows2+ctr3+2,3));
  fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
      ,'RS','H',ctr4+3,xx(nrows2+ctr3+3,1),...
      xx(nrows2+ctr3+3,2),xx(nrows2+ctr3+3,3));
  ctr3=ctr3+3;
  ctr4=ctr4+3;
    
end
ctr3=0;
%for ra=1:(nrows2/(natoms+1))
    
%  fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ra...
 %     ,strcat('RS',num2str(ra)),'H',nrows+ctr3+1,xx(nrows2+ctr3+1,1),...
 %     xx(nrows2+ctr3+1,2),xx(nrows2+ctr3+1,3));
  
 % fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ra...
 %     ,strcat('RS',num2str(ra)),'H',nrows+ctr3+2,xx(nrows2+ctr3+2,1),...
 %     xx(nrows2+ctr3+2,2),xx(nrows2+ctr3+2,3));
  
%  fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ra...
%      ,strcat('RS',num2str(ra)),'H',nrows+ctr3+3,xx(nrows2+ctr3+3,1),...
%      xx(nrows2+ctr3+3,2),xx(nrows2+ctr3+3,3));
  
%  ctr3=ctr3+3;
%end
fprintf(fid,'%7.5f %7.5f %7.5f',6.92,6.84,10);
fclose(fid);
xx;
