%CT - Any sp3 hybridised Carbon
%HC - Any aliphatic Hydrogen
%CT-HC=0.10900
leng = 8;
bread = 10;
lconst = 0.4943;
natoms = 11;  %Number of C-atoms
nrows1 = leng*bread; 
nrows2 = (natoms+1)*nrows1; %Number of C-atoms + Sulfur
nrows = (natoms+4)*nrows1; %Total number of atoms
xx=zeros(nrows,3);
offsetz = 3.5;
p = 1;

%Odd layers
for i=0:sqrt(3):(bread/2-1)*sqrt(3)
    for j = 1:leng
        xx(p,:) = [j*lconst i*lconst offsetz];  %Sulfur
        cr1 = 0;
        for k = 1:2:(natoms)
            xx(p+k,:) = [(j-0.08387)*lconst (i+cr1*0.0646)*lconst...
                offsetz+k*0.12916];
            cr1=cr1+1;
            
          %  xx(p+k+1,:)=[(j-0.08387)*lconst (i+cr1*0.0646)*lconst...
          %     -offsetz+k*0.12916];
            
        end
        cr2=0;
        for (l=2:2:natoms)
            xx(p+l,:)=[j*lconst (i+cr2*0.0646)*lconst... 
                offsetz+l*0.12916];
            
          % xx(p+l+1,:)=[j*lconst (i+cr2*0.0646)*lconst... 
          %     -offsetz+l*0.12916];
            
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
            
            
         %   xx(p+m+1,:)=[(j+0.5-0.08387)*lconst (i+cr3*0.0646)*lconst... 
         %       -offsetz+m*0.12916];
            cr3=cr3+1;
        end
        cr4=0;
        for (n=2:2:natoms)
            xx(p+n,:)=[(j+0.5)*lconst (i+cr4*0.0646)*lconst... 
                offsetz+n*0.12916];
                      
          %  xx(p+n+1,:)=[(j+0.5)*lconst (i+cr4*0.0646)*lconst... 
          %      offsetz+n*0.12916];
            cr4=cr4+1;
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


apratim = strcat('sam_methyl_both','.gro');
fid = fopen(apratim,'W');
fprintf(fid,'%s\n', 'SAM using Methyl head');
fprintf(fid,'%d\n',nrows1*29);
ctr=0;
ctr3=0;
ctr4=0; %Counter for atom number
translate = 0;
for ap=1:(nrows2/(natoms+1))
    ctr = ctr + 1;
    ctr4 = ctr4 + 1;
    fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
      ,'RS1','SC',ctr4,xx(ctr,1),xx(ctr,2),xx(ctr,3));
  ctr5=0;
    for c = ap + 1:ap+natoms
        ctr = ctr + 1;
        ctr5 = ctr5 + 1;
       
        fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
        ,'RS1',strcat('C',num2str(ctr5)),ctr4+1,xx(ctr,1),xx(ctr,2),xx(ctr,3));
              
        
       ctr4 = ctr4 + 1;
        
    end
  fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
      ,'RS1','Hml',ctr4+1,xx(nrows2+ctr3+1,1),...
      xx(nrows2+ctr3+1,2),xx(nrows2+ctr3+1,3));
  fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
      ,'RS1','Hml',ctr4+2,xx(nrows2+ctr3+2,1),...
      xx(nrows2+ctr3+2,2),xx(nrows2+ctr3+2,3));
  fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
      ,'RS1','Hml',ctr4+3,xx(nrows2+ctr3+3,1),...
      xx(nrows2+ctr3+3,2),xx(nrows2+ctr3+3,3));
  ctr4=ctr4+3;  
    ctr=ctr-11;
    
    ctr5=0;
    for c = ap + 1:ap+natoms
        ctr = ctr + 1;
        ctr5 = ctr5 + 1;
        ctr4=ctr4+1;
      fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
        ,'RS1',strcat('C',num2str(ctr5)),ctr4,xx(ctr,1),-xx(ctr,2)+translate*sqrt(3)*lconst...
        ,-xx(ctr,3)+2*offsetz);
        
    end
    
    
        
  
  fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
      ,'RS1','Hml',ctr4+1,xx(nrows2+ctr3+1,1),...
      -xx(nrows2+ctr3+1,2)+translate*sqrt(3)*lconst,-xx(nrows2+ctr3+1,3)+2*offsetz);
  fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
      ,'RS1','Hml',ctr4+2,xx(nrows2+ctr3+2,1),...
      -xx(nrows2+ctr3+2,2)+translate*sqrt(3)*lconst,-xx(nrows2+ctr3+2,3)+2*offsetz);
  fprintf(fid,'%5d%5s%5s%5d%8.3f%8.3f%8.3f\n',ap...
      ,'RS1','Hml',ctr4+3,xx(nrows2+ctr3+3,1),...
      -xx(nrows2+ctr3+3,2)+translate*sqrt(3)*lconst,-xx(nrows2+ctr3+3,3)+2*offsetz);
  
  ctr3=ctr3+3;
  ctr4=ctr4+3;
  if (mod(ap,leng)==0)
    translate=translate+1;
  end
    
end
ctr3=0;
fprintf(fid,'%7.5f %7.5f %7.5f',3.9544,4.2800,7);
fclose(fid);
xx;