function x = Hamins_ReadCH4(infile,outfile,mintime,ncol,nR,nV,R,Z)

cd(fileparts(mfilename('fullpath')))
cd ..
cd ..
cd ..
cd Validation/Hamins_CH4/FDS_Output_Files

c1 = 2;
c2 = c1 + nR;
c3 = c2 + nR;
c4 = c3 + nV;
c5 = c4 + nV;

RaR = zeros(1,nR);
RaV = zeros(1,nV);

if exist(fullfile(cd,infile), 'file') == 2
    fileID = fopen(fullfile(cd,infile));
    Reading = 'Reading';
  % display([Reading ' ' infile]);
    fgets(fileID);
    fgets(fileID);
    nt = 0;
    
    while(feof(fileID)==0)
        temp = str2num(fgets(fileID));
        M = zeros(1,ncol);
        M = temp;
        
        if (M(1)>=mintime)
            nt = nt + 1;
            Ra1 = [];
            Ra2 = [];
            Ra1(1:nR) = M(c1:c2-1);
            Ra2(1:nR) = M(c2:c3-1);
                        
            RaR = RaR + 0.5*(Ra1+Ra2);
            Ra1 = [];
            Ra2 = [];
            Ra1(1:nV) = M(c3:c4-1);
            Ra2(1:nV) = M(c4:c5-1);
            
            RaV = RaV + 0.5*(Ra1+Ra2);
        end
    end
    fclose(fileID);
    
    RaR = RaR / nt;
    RaV = RaV / nt;
  
else
     missing = 'is missing!';
     display([infile ' ' missing]);
     x=0;
     cd(fileparts(mfilename('fullpath')))
     return
end

titleString = 'Radius,Radial Flux,Height,Vertical Flux';
fileID = fopen(outfile,'w+');
fprintf(fileID, '%s', titleString);

for i = 1:nV
    if i <= nR
        fprintf(fileID, '\n %10.3e,%10.3e,%10.3e,%10.3e',R(i),RaR(i),Z(i),RaV(i));
    else
        fprintf(fileID, '\n %s,%10.3e,%10.3e','NaN,NaN',Z(i),RaV(i));
    end
end

fclose(fileID);
x=1;
cd(fileparts(mfilename('fullpath')))
return;
