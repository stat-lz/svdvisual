function svd3dplot(data, paramstruct);
% SVD3DPLOT displays Surface plot, image plot, SVD curve movie 
%                                       for Singular Value Decomposition
%
%         2004 Copyright (C) Lingsong Zhang (LSZHANG@email.unc.edu)
%
% Usage:
%          svd3dplot(data, parastruct);
%
% Inputs:
%
%   data         should be a MxN matrix, otherwise, no need to generate 3D
%                plot
%
%   paramstruct:
%         
%      a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1,...
%                            'field2',values2,...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%    fields             Value
%
%    icol               1 subtract the column mean first then do svd
%                         decomposition
%                       0 do not subtract the column mean(default)
%
%    irow               1 subtract the row mean first then do svd
%                         decomposition
%                       0 do not substract the row mean(default
%
%                       Note: 
%                       suggest only subtract one mean, either row mean or
%                       column mean. if both icol and irow are 1, we will
%                       subtract the double mean , i.e., the sum of the row
%                       mean and column mean, but  subtract the overall
%                       mean
%
%    idouble            0  do not subtract the double mean (default)
%
%                       1  subtract the double mean, the same as the icol
%                          as 1 and irow as 1. If idouble is given as 1, we
%                          will cover the icol and irow to be both 1. if
%                          idouble is given as 0, we will keep the same
%                          option of icol and irow. i.e. if you give icol
%                          as 1, irow as 1 and idouble as 0, we will still
%                          get a double mean subtraction.%                      
%
%    ioverall           0  do not subtract the overall mean (default)
%
%                       1  subtract the overall mean. If any one of the
%                          options irow, icol or idouble is 1, the program
%                          will discard this option.
%  
%    isamescale         0 Plot all the plots with tight scale, i.e.
%                         different plot using different scale
%
%                       1 Plot all the plots under the same scale. 
%
%                      
%
%    nosvd              number of SVD components (default is 3, if no mean
%                       subtraction; and 2, if any mean subtraction)
%
%    ips                1 generate ps output
%                       0 no ps output 
%
%    imovie             1 generate movie for svd components and double mean
%                         if there is.
%                       0 no movie generate (default)
%
%    iimage              0 no image view of the surface plots (default)
%
%                       1 image view of the surface plots
%
%    savestr            the file name you want to save. (default is svd surface plots)
%
%    axessize           The fontsize of axes labels and ticks, default is
%                       16
%
%    titlesize          the fontsize of subplot tiles
%
%    xlabelstr          String for x axis label (column direction), default
%                       is empty
%
%    ylabelstr          String for y axis label (row direction), default is
%                       empty
%
%
%    icolormap          0 (default) use the default color map in the matlab
%                          system. In this setting, red denotes maximum
%                          value and blue denotes minimum value
%
%                       1  use the gray color map, which generate plots
%                          for paper publication. In this setting, black
%                          denotes the minimum value and white denotes
%                          maximum value
%
%                       2  use the inverse-gray color map, which generate
%                          plots for paper publication. In this setting,
%                          black denotes the maximum value and white
%                          denotes minimum value.
%
% Outputs:
%     surface plots display (or ps file, if you give savestr input)
%     
%     no substract mean.  1-origianl data mesh plot
%                         2-first several svd components mesh plot
%                         3-first several svd components summation mesh plot
%                         4-residual mesh plot
%     number of svd components is no greater than 6.
%
%     with subtract mean.  1-origianl data mesh plot
%                          2-mean (column, row, or double) mesh plot 
%                          3-first several svd components mesh plot
%                          4-first several svd components summation mesh plot
%                          5-residual mesh plot
%     number of svd components is no greater than 5.
%
%
%Reminder: Please make sure svds.m in MatLab can be used.
%      Make sure the following files are in the same directory
%
%      svdls.m, columnmean.m, rowmean.m, doublemean.m, imagels.m by Lingsong Zhang
%
%(c)Copyright 2004, Lingsong Zhang (UNC-CH, STATISTICS DEPARTMENT)
%  revised by Lingsong Zhang (c) 2006
%Email: LSZHANG@email.unc.edu


icol=0;
irow=0;
ips=0;
nosvd=3;
imovie=0;
savestr='svd surface plots';
fs=16;
titlesize=11;
isamescale=0;
iimage=0;
idouble=0;
ioverall=0;
xlabelstr=' ';
ylabelstr=' ';
labelsize=11;
icolormap=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following block is used to input parameters                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin>1; %input parameters
  if isfield(paramstruct,'icol') ;    %  then change to input value
    icol = getfield(paramstruct,'icol') ; 
    nosvd = 2;
  end ; %end updating icol

  if isfield(paramstruct,'irow') ;    %  then change to input value
    irow = getfield(paramstruct,'irow') ; 
    nosvd = 2;
  end ; %end updating irow
  
  if isfield(paramstruct, 'nosvd');
      nosvd = getfield(paramstruct, 'nosvd');
  end;

  if isfield(paramstruct,'ips') ;    %  then change to input value
    ips = getfield(paramstruct,'ips') ; 
  end ; %end updating ips
  
  if isfield(paramstruct,'savestr') ;    %  then change to input value
    savestr = getfield(paramstruct,'savestr') ; 
  end ; %end updating savestr
  
  if isfield(paramstruct, 'imovie') ;
    imovie = getfield(paramstruct, 'imovie');
  end;
  
  if isfield(paramstruct, 'axessize');
     fs = getfield(paramstruct, 'axessize');
  end;
  
  if isfield(paramstruct, 'titlesize');
      titlesize = getfield(paramstruct, 'titlesize');
  end;
  
  if isfield(paramstruct, 'isamescale');
      isamescale = getfield(paramstruct, 'isamescale');
  end;
  
  if isfield(paramstruct, 'iimage');
      iimage = getfield(paramstruct, 'iimage');
  end;
  
  if isfield(paramstruct, 'idouble');
      idouble = getfield(paramstruct, 'idouble');
  end;
  
  if isfield(paramstruct, 'xlabelstr');
      xlabelstr = getfield(paramstruct, 'xlabelstr');
  end;
  
  if isfield(paramstruct, 'ylabelstr');
      ylabelstr =  getfield(paramstruct, 'ylabelstr');
  end;
  
  if isfield(paramstruct, 'labelsize');
      labelsize = getfield(paramstruct, 'labelsize');
  end;
  
  if isfield(paramstruct, 'icolormap');
      icolormap=getfield(paramstruct, 'icolormap');
  end;
  
  if isfield(paramstruct, 'ioverall');
      ioverall=getfield(paramstruct, 'ioverall');
  end;
end; %end inputting parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The above block is used to input parameters                    %
%                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[nrow, ncol]=size(data); %return the dimensions%
datamax=max(max(data));
datamin=min(min(data));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following block is used to design a special colormap       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colorvector=(1:-.01:0)';
mycolormap=[colorvector colorvector colorvector];
if nrow==1 | ncol==1;
    error('NO NEED TO GENERATE 3D PLOT!!!');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following block is used to check options irow, icol, idouble or ioverall   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if icol==1;
    ioverall=0;
end;%change the column mean option

if irow==1;
    ioverall=0;
end;%change the row mean option


if idouble==1;
    irow=1;
    icol=1;
    ioverall=0;
end;%change the double mean option.

if ioverall==1;
    irow=2;
    icol=2;
    %this is to overrule without removing any mean
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following block is used to generate the surface plot, image plot and more       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nosvd==0; %no svd component is needed.
  if irow==0 & icol==0;
      error('RAW SVD has no component!!! Check the option');

  elseif irow==0 & icol==1; %column svd
      colmat=columnmean(data);
      res=data-colmat;
     datass=sum(sum(data.^2));
     cmss=sum(sum(colmat.^2));
     cmp=cmss./datass;
     appp=cmp;
     
     figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;

      subplot(2, 2, 1);
          if iimage==0;
              mesh(data);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(data);
              set(gca, 'FontSize', fs);
          end;    
          title('Original data', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

      subplot(2, 2, 2);
          if iimage==0;
              mesh(colmat);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(colmat);
              set(gca, 'FontSize', fs);
          end;    
          title(['Column mean, ' num2str(cmp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

      subplot(2, 2, 3);
          if iimage==0;
              mesh(colmat);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(colmat);
              set(gca, 'FontSize', fs);
          end;    
          title('Approximation', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

      subplot(2, 2, 4);
          if iimage==0;
              mesh(res);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(res);
              set(gca, 'FontSize', fs);
          end;    
          title('Residual', 'FontSize', titlesize);          
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
      
  elseif irow==1 & icol==0;
      rowmat=rowmean(data);
      res=data-rowmat;
     datass=sum(sum(data.^2));
     rmss=sum(sum(rowmat.^2));
     rmp=rmss./datass;
     appp=rmp;

     figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;
     
     
      subplot(2, 2, 1);
          if iimage==0;
              mesh(data);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(data);
              set(gca, 'FontSize', fs);
          end;    
          title('Original data', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

      subplot(2, 2, 2);
          if iimage==0;
              mesh(rowmat);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(rowmat);
              set(gca, 'FontSize', fs);
          end;    
          title(['Row mean, ' num2str(rmp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

      subplot(2, 2, 3);
          if iimage==0;
              mesh(rowmat);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(rowmat);
              set(gca, 'FontSize', fs);
          end;    
          title('Approximation', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

      subplot(2, 2, 4);
          if iimage==0;
              mesh(res);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(res);
              set(gca, 'FontSize', fs);
          end;    
          title('Residual', 'FontSize', titlesize);          
           xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
           
      
  elseif irow==1 & icol==1;
           doubmat=doublemean(data);
      res=data-doubmat;
     datass=sum(sum(data.^2));
     dmss=sum(sum(doubmat.^2));
     dmp=dmss./datass;
     appp=dmp;

     figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;

     subplot(2, 2, 1);
          if iimage==0;
              mesh(data);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(data);
              set(gca, 'FontSize', fs);
          end;    
          title('Original data', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

          subplot(2, 2, 2);
          if iimage==0;
              mesh(doubmat);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(doubmat);
              set(gca, 'FontSize', fs);
          end;    
          title(['Double mean, ' num2str(dmp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

      subplot(2, 2, 3);
          if iimage==0;
              mesh(doubmat);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(doubmat);
              set(gca, 'FontSize', fs);
          end;    
          title('Approximation', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

      subplot(2, 2, 4);
          if iimage==0;
              mesh(res);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(res);
              set(gca, 'FontSize', fs);
          end;    
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
          title('Residual', 'FontSize', titlesize); 
%  end;
  
%overall mean

  elseif ioverall==1;
     omeanmat=overallmeanmean(data);
     res=data-omeanmat;
     datass=sum(sum(data.^2));
     omss=sum(sum(omeanmat.^2));
     omp=omss./datass;
     appp=omp;

     figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;

     subplot(2, 2, 1);
          if iimage==0;
              mesh(data);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(data);
              set(gca, 'FontSize', fs);
          end;    
          title('Original data', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

          subplot(2, 2, 2);
          if iimage==0;
              mesh(omeanmat);
              set(gca, 'FontSize', fs);
%               if isamescale==0;
%                   axis tight;
%               elseif isamescale==1;
              axis([0 ncol 0 nrow datamin datamax]);
%               end;
          else
              imagels(omeanmat);
              set(gca, 'FontSize', fs);
          end;    
          title(['Overall mean, ' num2str(omp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

      subplot(2, 2, 3);
          if iimage==0;
              mesh(omeanmat);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(omeanmat);
              set(gca, 'FontSize', fs);
          end;    
          title('Approximation', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

      subplot(2, 2, 4);
          if iimage==0;
              mesh(res);
              set(gca, 'FontSize', fs);
              if isamescale==0;
                  axis tight;
              elseif isamescale==1;
                  axis([0 ncol 0 nrow datamin datamax]);
              end;
          else
              imagels(res);
              set(gca, 'FontSize', fs);
          end;    
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
          title('Residual', 'FontSize', titlesize); 
  end;

%end overall mean
  
else %at least one svd component is needed

  if irow==0 & icol==0;
    if nosvd>7;
        error('NUMBER OF SVD COMPONENTS CAN NOT BE GREATER THAN 6!!!');
    end;

    disp('No mean subtraction!');%indicate mean subtraction
    
    [u, s, v]=svdls(data, nosvd);
    app=u*s*v';
    res=data-app; %finish calculation
    
     figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;
    
    if nosvd==1;
        nprow=2;
        npcol=2;
        nbegin=1;
        napp=3;
        nres=4;
    elseif nosvd==2;
        nprow=3;
        npcol=2;
        nbegin=2;
        napp=5;
        nres=6;
    elseif nosvd==3;
        nprow=3;
        npcol=2;
        nbegin=1;
        napp=5;
        nres=6;
    elseif nosvd==4;
        nprow=3;
        npcol=3;
        nbegin=2;
        napp=7;
        nres=9;
    elseif nosvd==5;
        nprow=3;
        npcol=3;
        nbegin=1;
        napp=7;
        nres=9;
    else
        nprow=3;
        npcol=3;
        nbegin=1;
        napp=8;
        nres=9;
    end;

    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;

    clf;%display original data matrix
    subplot(nprow, npcol, 1);
    if iimage==0;
    mesh(data);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;    
    else
        imagels(data);
        set(gca, 'FontSize', fs);
    end;
    title('Original data', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
        
    clear datass;
    datass=sum(sum(data.^2));
    appss=sum(sum(app.^2));
    appp=appss./datass;
    
    for i=1:nosvd; %display svd components
        SVtemp=s(i, i)*u(:, i)*v(:, i)';
        SVtempSS=sum(sum(SVtemp.^2));
        SVtempp=SVtempSS./datass;
        subplot(nprow, npcol, nbegin+i);
        if iimage==0;
        mesh(SVtemp);
        set(gca, 'FontSize', fs);
        if isamescale==0;
            axis tight;
        elseif isamescale==1;
            axis([0 ncol 0 nrow datamin datamax]);
        end;
        else
            imagels(SVtemp);
            set(gca, 'FontSize', fs);
        end;

        title(['SV' num2str(i) ', ' num2str(SVtempp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
    end;
    
    subplot(nprow, npcol, napp);%display the first several approximation
    if iimage==0;
    mesh(app);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
    else
        imagels(app);
        set(gca, 'FontSize', fs);
    end;
    title(['First ' num2str(nosvd) ' SVs approximation, ' num2str(appp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
    
    subplot(nprow, npcol, nres);%display the residual component
    if iimage==0;
    mesh(res);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(res);
        set(gca, 'FontSize', fs);
    end;


    title('Residual', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
    
elseif irow==0 & icol==1;
    
     figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;

    if nosvd>5;
        error('NUMBER OF SVD COMPONENTS CAN NOT BE GREATER THAN 5!!!');
    end;

    disp('Column mean subtraction!');%indicate mean subtraction
    colmat=columnmean(data);
    locdata=data-colmat;
    
    [u, s, v]=svdls(locdata, nosvd);
    svdapp=u*s*v';
    app=colmat+svdapp;
    res=data-app; %finish calculation
    
    
    if nosvd==1;
        nprow=3;
        npcol=2;
        nbegin=2;
        napp=5;
        nres=6;
    elseif nosvd==2;
        nprow=3;
        npcol=2;
        nbegin=2;
        napp=5;
        nres=6;
    elseif nosvd==3;
        nprow=3;
        npcol=3;
        nbegin=3;
        napp=7;
        nres=9;
    elseif nosvd==4;
        nprow=3;
        npcol=3;
        nbegin=3;
        napp=8;
        nres=9;
    elseif nosvd==5;
        nprow=3;
        npcol=3;
        nbegin=2;
        napp=8;
        nres=9;
    end;

    clf;
    subplot(nprow, npcol, 1);
    if iimage==0;
    mesh(data);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
    else
        imagels(data);
        set(gca, 'FontSize', fs);
    end;


    title('Original data', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

    clear datass;
    datass=sum(sum(data.^2));
    appss=sum(sum(app.^2));
    cmss=sum(sum(colmat.^2));
    cmp=cmss./datass;
    appp=appss./datass;
    
    subplot(nprow, npcol, 2);
    if iimage==0;
    mesh(colmat);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
    else
        imagels(colmat);
        set(gca, 'FontSize', fs);
    end;


    title(['Column mean, ' num2str(cmp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

    for i=1:nosvd;
        SVtemp=s(i, i)*u(:, i)*v(:, i)';
        SVtempSS=sum(sum(SVtemp.^2));
        SVtempp=SVtempSS./datass;
        subplot(nprow, npcol, nbegin+i);
        if iimage==0;
        mesh(SVtemp);
        set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(SVtemp);
        set(gca, 'FontSize', fs);
    end;


          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

    title(['SV' num2str(i) ', ' num2str(SVtempp) ' of TSS'], 'FontSize', titlesize);
    end;
    
    subplot(nprow, npcol, napp);
    if iimage==0;
    mesh(app);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
    else
        imagels(app);
        set(gca, 'FontSize', fs);
    end;

    title(['Approximation, ' num2str(appp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

         subplot(nprow, npcol, nres);
    if iimage==0;
    mesh(res);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
    else
        imagels(res);
        set(gca, 'FontSize', fs);
    end;

    title('Residual', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
    

elseif irow==1 & icol==0;
    if nosvd>5;
        error('NUMBER OF SVD COMPONENTS CAN NOT BE GREATER THAN 5!!!');
    end;

    disp('Row mean subtraction!');%indicate mean subtraction
    rowmat=rowmean(data);
    locdata=data-rowmat;
    
    [u, s, v]=svdls(locdata, nosvd);
    svdapp=u*s*v';
    app=rowmat+svdapp;
    res=data-app; %finish calculation
    
     figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;
    
    if nosvd==1;
        nprow=3;
        npcol=2;
        nbegin=2;
        napp=5;
        nres=6;
    elseif nosvd==2;
        nprow=3;
        npcol=2;
        nbegin=2;
        napp=5;
        nres=6;
    elseif nosvd==3;
        nprow=3;
        npcol=3;
        nbegin=3;
        napp=7;
        nres=9;
    elseif nosvd==4;
        nprow=3;
        npcol=3;
        nbegin=3;
        napp=8;
        nres=9;
    elseif nosvd==5;
        nprow=3;
        npcol=3;
        nbegin=2;
        napp=8;
        nres=9;
    end;

    clf;
    subplot(nprow, npcol, 1);
    if iimage==0;
    mesh(data);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
    else
        imagels(data);
        set(gca, 'FontSize', fs);
    end;


    title('Original data', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

    clear datass;
    datass=sum(sum(data.^2));
    appss=sum(sum(app.^2));
    rmss=sum(sum(rowmat.^2));
    rmp=rmss./datass;
    appp=appss./datass;
    
    subplot(nprow, npcol, 2);
    if iimage==0;
    mesh(rowmat);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
    else
        imagels(rowmat);
        set(gca, 'FontSize', fs);
    end;


    title(['Row mean, ' num2str(rmp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

    for i=1:nosvd;
        SVtemp=s(i, i)*u(:, i)*v(:, i)';
        SVtempSS=sum(sum(SVtemp.^2));
        SVtempp=SVtempSS./datass;
        subplot(nprow, npcol, nbegin+i);
        if iimage==0;
        mesh(SVtemp);
        set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(SVtemp);
        set(gca, 'FontSize', fs);
    end;


        title(['SV' num2str(i) ', ' num2str(SVtempp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

    end;
    
    subplot(nprow, npcol, napp);
    if iimage==0;
    mesh(app);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(app);
        set(gca, 'FontSize', fs);
    end;


    title(['Approximation, ' num2str(appp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

          subplot(nprow, npcol, nres);
    if iimage==0;
    mesh(res);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(res);
        set(gca, 'FontSize', fs);
    end;

    title('Residual', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
    

elseif irow==1 & icol==1;
    if nosvd>5;
        error('NUMBER OF SVD COMPONENTS CAN NOT BE GREATER THAN 5!!!');
    end;

    disp('Double mean subtraction!');%indicate mean subtraction
    doubmat=doublemean(data);
    locdata=data-doubmat;
    
    [u, s, v]=svdls(locdata, nosvd);
    svdapp=u*s*v';
    app=doubmat+svdapp;
    res=data-app; %finish calculation
    
     figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;
    
    if nosvd==1;
        nprow=3;
        npcol=2;
        nbegin=2;
        napp=5;
        nres=6;
    elseif nosvd==2;
        nprow=3;
        npcol=2;
        nbegin=2;
        napp=5;
        nres=6;
    elseif nosvd==3;
        nprow=3;
        npcol=3;
        nbegin=3;
        napp=7;
        nres=9;
    elseif nosvd==4;
        nprow=3;
        npcol=3;
        nbegin=3;
        napp=8;
        nres=9;
    elseif nosvd==5;
        nprow=3;
        npcol=3;
        nbegin=2;
        napp=8;
        nres=9;
    end;

    clf;
    subplot(nprow, npcol, 1);
    if iimage==0;
    mesh(data);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(data);
        set(gca, 'FontSize', fs);
    end;


    title('Original data', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
    
    clear datass;
    datass=sum(sum(data.^2));
    appss=sum(sum(app.^2));
    dmss=sum(sum(doubmat.^2));
    dmp=dmss./datass;
    appp=appss./datass;
    
    subplot(nprow, npcol, 2);
    if iimage==0;
    mesh(doubmat);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(doubmat);
        set(gca, 'FontSize', fs);
    end;

    title(['Double mean, ' num2str(dmp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

    for i=1:nosvd;
        SVtemp=s(i, i)*u(:, i)*v(:, i)';
        SVtempSS=sum(sum(SVtemp.^2));
        SVtempp=SVtempSS./datass;
        subplot(nprow, npcol, nbegin+i);
        if iimage==0;
        mesh(SVtemp);
        set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(SVtemp);
        set(gca, 'FontSize', fs);
    end;

          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
        title(['SV' num2str(i) ', ' num2str(SVtempp) ' of TSS'], 'FontSize', titlesize);
    end;
    
    subplot(nprow, npcol, napp);
    if iimage==0;
    mesh(app);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(app);
        set(gca, 'FontSize', fs);
    end;

    title(['Approximation, ' num2str(appp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

          subplot(nprow, npcol, nres);
    if iimage==0;
    mesh(res);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(res);
        set(gca, 'FontSize', fs);
    end;

    title('Residual', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

%  end;

  %overall mean should be removed
  elseif ioverall==1;
    if nosvd>5;
        error('NUMBER OF SVD COMPONENTS CAN NOT BE GREATER THAN 5!!!');
    end;

    disp('Overall mean subtraction!');%indicate mean subtraction
    omeanmat=overallmean(data);
    locdata=data-omeanmat;
    
    [u, s, v]=svdls(locdata, nosvd);
    svdapp=u*s*v';
    app=omeanmat+svdapp;
    res=data-app; %finish calculation
    
     figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;
    
    if nosvd==1;
        nprow=3;
        npcol=2;
        nbegin=2;
        napp=5;
        nres=6;
    elseif nosvd==2;
        nprow=3;
        npcol=2;
        nbegin=2;
        napp=5;
        nres=6;
    elseif nosvd==3;
        nprow=3;
        npcol=3;
        nbegin=3;
        napp=7;
        nres=9;
    elseif nosvd==4;
        nprow=3;
        npcol=3;
        nbegin=3;
        napp=8;
        nres=9;
    elseif nosvd==5;
        nprow=3;
        npcol=3;
        nbegin=2;
        napp=8;
        nres=9;
    end;

    clf;
    subplot(nprow, npcol, 1);
    if iimage==0;
    mesh(data);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(data);
        set(gca, 'FontSize', fs);
    end;


    title('Original data', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
    
    clear datass;
    datass=sum(sum(data.^2));
    appss=sum(sum(app.^2));
    omss=sum(sum(omeanmat.^2));
    omp=omss./datass;
    appp=appss./datass;
    
    subplot(nprow, npcol, 2);
    if iimage==0;
    mesh(omeanmat);
    set(gca, 'FontSize', fs);
%     if isamescale==0;
%         axis tight;
%     elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
%     end;
        else
        imagels(omeanmat);
        set(gca, 'FontSize', fs);
    end;

    title(['Overall mean, ' num2str(omp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

    for i=1:nosvd;
        SVtemp=s(i, i)*u(:, i)*v(:, i)';
        SVtempSS=sum(sum(SVtemp.^2));
        SVtempp=SVtempSS./datass;
        subplot(nprow, npcol, nbegin+i);
        if iimage==0;
        mesh(SVtemp);
        set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(SVtemp);
        set(gca, 'FontSize', fs);
    end;

          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
        title(['SV' num2str(i) ', ' num2str(SVtempp) ' of TSS'], 'FontSize', titlesize);
    end;
    
    subplot(nprow, npcol, napp);
    if iimage==0;
    mesh(app);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(app);
        set(gca, 'FontSize', fs);
    end;

    title(['Approximation, ' num2str(appp) ' of TSS'], 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

          subplot(nprow, npcol, nres);
    if iimage==0;
    mesh(res);
    set(gca, 'FontSize', fs);
    if isamescale==0;
        axis tight;
    elseif isamescale==1;
        axis([0 ncol 0 nrow datamin datamax]);
    end;
        else
        imagels(res);
        set(gca, 'FontSize', fs);
    end;

    title('Residual', 'FontSize', titlesize);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);

  end;
  
  
  
  
end; %nosvd>0

if ips==1;
    orient tall;
    print('-dpsc', [savestr '.ps']);
end;


%the following is used to generate movies
if imovie==1;
    if nrow>100 | ncol>100;
        error('The dimensionality of data matrix is too large, zoomed movie is recommended!')
    end;

    if irow==0 & icol==0;
        disp('Movie generated for no mean subtraction');
        
        %set parameters
        moviefps = 5;
        moviecstr = 'Cinepak' ;
        nrepeat = 1;
        nframe = 1+nosvd+nosvd*(ncol+nrow);
        
        %initiate figure and movie
        fighand = figure(1);
        clf;
        clear moviestruct;
        
        %Do main movie loop
        %original data mesh plot;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;

     mesh(data);
        set(gca, 'FontSize', fs);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
        zlabel('Original data');
        moviestruct(1) = getframe(fighand);
        
        
        for isvd=1:nosvd;
            SVtemp=s(isvd, isvd)*u(:, isvd)*v(:, isvd)';
            Utemp=u(:, isvd);
            Vtemp=v(:, isvd);
            
            if Vtemp<0;
                Vtemp=-1.*Vtemp;
                Utemp=-1.*Utemp;
            end;
            
            Zmax=max(max(SVtemp));
            Zmin=min(min(SVtemp));
            Vshape=Zmax*Vtemp./max(abs(Vtemp));
            Ushape=Zmax*Utemp./max(abs(Utemp));
            
            %mesh plot for the SVD components
            colormap('default');
            mesh(SVtemp);
            set(gca, 'FontSize', fs);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
            zlabel(['SV' num2str(isvd)]);
            moviestruct((1+nrow+ncol)*(isvd-1)+2)=getframe(fighand);
            
            %set initial values
            x0=(1:ncol);
            y0=(1:nrow);
            yx0=zeros(1, ncol);
            xy0=zeros(1, nrow);
            yx=yx0;
            xy=xy0;
            
            for irow=1:nrow;
                yx=yx+1;
                
                colormap('gray');
                mesh(SVtemp);
                set(gca, 'FontSize', fs);
                xlabel(xlabelstr, 'FontSize', labelsize);
                ylabel([ylabelstr, ' ', num2str(irow)], 'FontSize', labelsize);
                zlabel(['SV', num2str(isvd)]);
                hold on;
                plot3(x0, yx0, Vshape, 'b-.', 'LineWidth', 4);
                plot3(xy0, y0, Ushape, 'g--', 'LineWidth', 4);
                plot3(x0, yx, SVtemp(irow, :), 'r-', 'LineWidth', 6);
                plot3(xy0(irow), y0(irow), Ushape(irow), 'r.', 'MarkerSize', 50);
                hold off;
                
                moviestruct((1+nrow+ncol)*(isvd-1)+2+irow) = getframe(fighand);
            end;
            
            for icol=1:ncol;
                xy=xy+1;
                
                colormap('gray');
                mesh(SVtemp);
                set(gca, 'FontSize', fs);
                xlabel([xlabelstr, ' ', num2str(icol)], 'FontSize', labelsize);
                ylabel(ylabelstr, 'FontSize', labelsize);
                zlabel(['SV', num2str(isvd)]);
                hold on;
                plot3(x0, yx0, Vshape, 'b-.', 'LineWidth', 4);
                plot3(xy0, y0, Ushape, 'g--', 'LineWidth', 4);
                plot3(xy, y0, SVtemp(:, icol), 'r-', 'LineWidth', 6);
                plot3(x0(icol), yx0(icol), Vshape(icol), 'r.', 'MarkerSize', 50);
                hold off;
                
                moviestruct((1+nrow+ncol)*(isvd-1)+2+nrow+icol) = getframe(fighand);
            end;
        end;

        fullorder=[(1:nrow), (nrow-1:-1:1), (1+nrow:nrow+ncol), ((nrow+ncol-1):-1:1+nrow)];
            
        initone=ones(1, 5);
        clear isvd;
        vmorder=initone;
        
        for isvd=1:nosvd;
            firstposition=(1+nrow+ncol)*(isvd-1)+2;
            innerone=firstposition*initone;
            innerposition=firstposition+fullorder;
            temp=vmorder;
            vmorder=[temp innerone innerposition];
        end;
      
        moviestruct = moviestruct(vmorder);
            
        savestr='SVDmovie';
            
        movie2avi(moviestruct,savestr,'compression',moviecstr, ...
                      'keyframe',moviefps,'fps',moviefps) ;

        
        
    
    
    end; %end irow and icol

    if irow==0 & icol==1;
        disp('Movie generated for taking column mean');
        
        %set parameters
        moviefps = 5;
        moviecstr = 'Cinepak' ;
        nrepeat = 1;
        nframe = 2+nosvd+nosvd*(ncol+nrow);
        
        %initiate figure and movie
        fighand = figure(1);
        clf;
        clear moviestruct;
        
        %Do main movie loop
        %original data mesh plot;
        
          figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;
     
     mesh(data);
        set(gca, 'FontSize', fs);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
        zlabel('Original data');
        moviestruct(1) = getframe(fighand);
        
        colormap('default');
        mesh(colmat);
        set(gca, 'FontSize', fs);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
        zlabel('Column mean');
        moviestruct(2) = getframe(fighand);
        
        
        for isvd=1:nosvd;
            SVtemp=s(isvd, isvd)*u(:, isvd)*v(:, isvd)';
            Utemp=u(:, isvd);
            Vtemp=v(:, isvd);
            
            if Vtemp<0;
                Vtemp=-1.*Vtemp;
                Utemp=-1.*Utemp;
            end;
            
            Zmax=max(max(SVtemp));
            Zmin=min(min(SVtemp));
            Vshape=Zmax*Vtemp./max(abs(Vtemp));
            Ushape=Zmax*Utemp./max(abs(Utemp));
            
            %mesh plot for the SVD components
            colormap('default');
            mesh(SVtemp);
            set(gca, 'FontSize', fs);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
            zlabel(['SV' num2str(isvd)]);
            moviestruct((1+nrow+ncol)*(isvd-1)+3)=getframe(fighand);
            
            %set initial values
            x0=(1:ncol);
            y0=(1:nrow);
            yx0=zeros(1, ncol);
            xy0=zeros(1, nrow);
            yx=yx0;
            xy=xy0;
            
            for irow=1:nrow;
                yx=yx+1;
                
                colormap('gray');
                mesh(SVtemp);
                set(gca, 'FontSize', fs);
                xlabel(xlabelstr, 'FontSize', labelsize);
                ylabel([ylabelstr, ' ', num2str(irow)], 'FontSize', labelsize);
                zlabel(['SV', num2str(isvd)]);
                hold on;
                plot3(x0, yx0, Vshape, 'b-.', 'LineWidth', 4);
                plot3(xy0, y0, Ushape, 'g--', 'LineWidth', 4);
                plot3(x0, yx, SVtemp(irow, :), 'r-', 'LineWidth', 6);
                plot3(xy0(irow), y0(irow), Ushape(irow), 'r.', 'MarkerSize', 50);
                hold off;
                
                moviestruct((1+nrow+ncol)*(isvd-1)+3+irow) = getframe(fighand);
            end;
            
            for icol=1:ncol;
                xy=xy+1;
                
                colormap('gray');
                mesh(SVtemp);
                set(gca, 'FontSize', fs);
                xlabel([xlabelstr, ' ', num2str(icol)], 'FontSize', labelsize);
                ylabel(ylabelstr, 'FontSize', labelsize);
                zlabel(['SV', num2str(isvd)]);
                hold on;
                plot3(x0, yx0, Vshape, 'b-.', 'LineWidth', 4);
                plot3(xy0, y0, Ushape, 'g--', 'LineWidth', 4);
                plot3(xy, y0, SVtemp(:, icol), 'r-', 'LineWidth', 6);
                plot3(x0(icol), yx0(icol), Vshape(icol), 'r.', 'MarkerSize', 50);
                hold off;
                
                moviestruct((1+nrow+ncol)*(isvd-1)+3+nrow+icol) = getframe(fighand);
            end;
        end;

        fullorder=[(1:nrow), (nrow-1:-1:1), (1+nrow:nrow+ncol), ((nrow+ncol-1):-1:1+nrow)];
            
        initone=ones(1, 5);
        clear isvd;
        vmorder=[initone, 2*initone];
        
        for isvd=1:nosvd;
            firstposition=(1+nrow+ncol)*(isvd-1)+3;
            innerone=firstposition*initone;
            innerposition=firstposition+fullorder;
            temp=vmorder;
            vmorder=[temp innerone innerposition];
        end;
      
        moviestruct = moviestruct(vmorder);
            
        savestr='SVDmovie';
            
        movie2avi(moviestruct,savestr,'compression',moviecstr, ...
                      'keyframe',moviefps,'fps',moviefps) ;
    end; %end irow and icol
    
    
    if irow==1 & icol==0;
        disp('Movie generated for taking row mean');
        
        %set parameters
        moviefps = 5;
        moviecstr = 'Cinepak' ;
        nrepeat = 1;
        nframe = 2+nosvd+nosvd*(ncol+nrow);
        
        %initiate figure and movie
        fighand = figure(1);
        clf;
        clear moviestruct;
        
        %Do main movie loop
        %original data mesh plot;
        
          figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;

     mesh(data);
        set(gca, 'FontSize', fs);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
        zlabel('Original data');
        moviestruct(1) = getframe(fighand);
        
        colormap('default');
        mesh(rowmat);
        set(gca, 'FontSize', fs);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
        zlabel('Row mean');
        moviestruct(2) = getframe(fighand);
        
        
        for isvd=1:nosvd;
            SVtemp=s(isvd, isvd)*u(:, isvd)*v(:, isvd)';
            Utemp=u(:, isvd);
            Vtemp=v(:, isvd);
            
            if Vtemp<0;
                Vtemp=-1.*Vtemp;
                Utemp=-1.*Utemp;
            end;
            
            Zmax=max(max(SVtemp));
            Zmin=min(min(SVtemp));
            Vshape=Zmax*Vtemp./max(abs(Vtemp));
            Ushape=Zmax*Utemp./max(abs(Utemp));
            
            %mesh plot for the SVD components
            colormap('default');
            mesh(SVtemp);
            set(gca, 'FontSize', fs);
          xlabel(xlabelstr, 'FontSize', labelsize);
          ylabel(ylabelstr, 'FontSize', labelsize);
            zlabel(['SV' num2str(isvd)]);
            moviestruct((1+nrow+ncol)*(isvd-1)+3)=getframe(fighand);
            
            %set initial values
            x0=(1:ncol);
            y0=(1:nrow);
            yx0=zeros(1, ncol);
            xy0=zeros(1, nrow);
            yx=yx0;
            xy=xy0;
            
            for irow=1:nrow;
                yx=yx+1;
                
                colormap('gray');
                mesh(SVtemp);
                set(gca, 'FontSize', fs);
                xlabel(xlabelstr, 'FontSize', labelsize);
                ylabel([ylabelstr, ' ', num2str(irow)], 'FontSize', labelsize);
                zlabel(['SV', num2str(isvd)]);
                hold on;
                plot3(x0, yx0, Vshape, 'b-.', 'LineWidth', 4);
                plot3(xy0, y0, Ushape, 'g--', 'LineWidth', 4);
                plot3(x0, yx, SVtemp(irow, :), 'r-', 'LineWidth', 6);
                plot3(xy0(irow), y0(irow), Ushape(irow), 'r.', 'MarkerSize', 50);
                hold off;
                
                moviestruct((1+nrow+ncol)*(isvd-1)+3+irow) = getframe(fighand);
            end;
            
            for icol=1:ncol;
                xy=xy+1;
                
                colormap('gray');
                mesh(SVtemp);
                set(gca, 'FontSize', fs);
                xlabel([xlabelstr, ' ', num2str(icol)], 'FontSize', labelsize);
                ylabel(ylabelstr, 'FontSize', labelsize);
                zlabel(['SV', num2str(isvd)]);
                hold on;
                plot3(x0, yx0, Vshape, 'b-.', 'LineWidth', 4);
                plot3(xy0, y0, Ushape, 'g--', 'LineWidth', 4);
                plot3(xy, y0, SVtemp(:, icol), 'r-', 'LineWidth', 6);
                plot3(x0(icol), yx0(icol), Vshape(icol), 'r.', 'MarkerSize', 50);
                hold off;
                
                moviestruct((1+nrow+ncol)*(isvd-1)+3+nrow+icol) = getframe(fighand);
            end;
        end;

        fullorder=[(1:nrow), (nrow-1:-1:1), (1+nrow:nrow+ncol), ((nrow+ncol-1):-1:1+nrow)];
            
        initone=ones(1, 5);
        clear isvd;
        vmorder=[initone, 2*initone];
        
        for isvd=1:nosvd;
            firstposition=(1+nrow+ncol)*(isvd-1)+3;
            innerone=firstposition*initone;
            innerposition=firstposition+fullorder;
            temp=vmorder;
            vmorder=[temp innerone innerposition];
        end;
      
        moviestruct = moviestruct(vmorder);
            
        savestr='SVDmovie';
            
        movie2avi(moviestruct,savestr,'compression',moviecstr, ...
                      'keyframe',moviefps,'fps',moviefps) ;

        
        
    
    
    end; %end irow and icol
    
    
    if irow==1 & icol==1;
        disp('Movie generated for taking double mean');
        
        %set parameters
        moviefps = 5;
        moviecstr = 'Cinepak' ;
        nrepeat = 1;
        nframe = 2+nosvd+nosvd*(ncol+nrow);
        
        %initiate figure and movie
        fighand = figure(1);
        clf;
        clear moviestruct;
        
        %Do main movie loop
        %original data mesh plot;
        
          figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;

     mesh(data);
        set(gca, 'FontSize', fs);
        xlabel(xlabelstr, 'FontSize', labelsize);
        ylabel(ylabelstr, 'FontSize', labelsize);
        zlabel('Original data');
        moviestruct(1) = getframe(fighand);
        
        colormap('default');
        mesh(doubmat);
        set(gca, 'FontSize', fs);
        xlabel(xlabelstr, 'FontSize', labelsize);
        ylabel(ylabelstr, 'FontSize', labelsize);
        zlabel('Double mean');
        moviestruct(2) = getframe(fighand);
        
        
        for isvd=1:nosvd;
            SVtemp=s(isvd, isvd)*u(:, isvd)*v(:, isvd)';
            Utemp=u(:, isvd);
            Vtemp=v(:, isvd);
            
            if Vtemp<0;
                Vtemp=-1.*Vtemp;
                Utemp=-1.*Utemp;
            end;
            
            Zmax=max(max(SVtemp));
            Zmin=min(min(SVtemp));
            Vshape=Zmax*Vtemp./max(abs(Vtemp));
            Ushape=Zmax*Utemp./max(abs(Utemp));
            
            %mesh plot for the SVD components
            colormap('default');
            mesh(SVtemp);
            set(gca, 'FontSize', fs);
        xlabel(xlabelstr, 'FontSize', labelsize);
        ylabel(ylabelstr, 'FontSize', labelsize);
            zlabel(['SV' num2str(isvd)]);
            moviestruct((1+nrow+ncol)*(isvd-1)+3)=getframe(fighand);
            
            %set initial values
            x0=(1:ncol);
            y0=(1:nrow);
            yx0=zeros(1, ncol);
            xy0=zeros(1, nrow);
            yx=yx0;
            xy=xy0;
            
            for irow=1:nrow;
                yx=yx+1;
                
                colormap('gray');
                mesh(SVtemp);
                set(gca, 'FontSize', fs);
                xlabel(xlabelstr, 'FontSize', labelsize);
                ylabel([ylabelstr, ' ', num2str(irow)], 'FontSize', labelsize);
                zlabel(['SV', num2str(isvd)]);
                hold on;
                plot3(x0, yx0, Vshape, 'b-.', 'LineWidth', 4);
                plot3(xy0, y0, Ushape, 'g--', 'LineWidth', 4);
                plot3(x0, yx, SVtemp(irow, :), 'r-', 'LineWidth', 6);
                plot3(xy0(irow), y0(irow), Ushape(irow), 'r.', 'MarkerSize', 50);
                hold off;
                
                moviestruct((1+nrow+ncol)*(isvd-1)+3+irow) = getframe(fighand);
            end;
            
            for icol=1:ncol;
                xy=xy+1;
                
                colormap('gray');
                mesh(SVtemp);
                set(gca, 'FontSize', fs);
                xlabel([xlabelstr, ' ', num2str(icol)], 'FontSize', labelsize);
                ylabel(ylabelstr, 'FontSize', labelsize);
                zlabel(['SV', num2str(isvd)]);
                hold on;
                plot3(x0, yx0, Vshape, 'b-.', 'LineWidth', 4);
                plot3(xy0, y0, Ushape, 'g--', 'LineWidth', 4);
                plot3(xy, y0, SVtemp(:, icol), 'r-', 'LineWidth', 6);
                plot3(x0(icol), yx0(icol), Vshape(icol), 'r.', 'MarkerSize', 50);
                hold off;
                
                moviestruct((1+nrow+ncol)*(isvd-1)+3+nrow+icol) = getframe(fighand);
            end;
        end;

        fullorder=[(1:nrow), (nrow-1:-1:1), (1+nrow:nrow+ncol), ((nrow+ncol-1):-1:1+nrow)];
            
        initone=ones(1, 5);
        clear isvd;
        vmorder=[initone, 2*initone];
        
        for isvd=1:nosvd;
            firstposition=(1+nrow+ncol)*(isvd-1)+3;
            innerone=firstposition*initone;
            innerposition=firstposition+fullorder;
            temp=vmorder;
            vmorder=[temp innerone innerposition];
        end;
      
        moviestruct = moviestruct(vmorder);
            
        savestr='SVDmovie';
            
        movie2avi(moviestruct,savestr,'compression',moviecstr, ...
                      'keyframe',moviefps,'fps',moviefps) ;

        
        
    
    
    end; %end irow and icol
   
        if ioverall==1;
        disp('Movie generated for taking overall mean');
        
        %set parameters
        moviefps = 5;
        moviecstr = 'Cinepak' ;
        nrepeat = 1;
        nframe = 2+nosvd+nosvd*(ncol+nrow);
        
        %initiate figure and movie
        fighand = figure(1);
        clf;
        clear moviestruct;
        
        %Do main movie loop
        %original data mesh plot;
        
          figure; clf;
     if icolormap==0;
         colormap('default');
     elseif icolormap==1;
         colormap('gray');
     elseif icolormap==2;
         colormap(mycolormap);
     end;

     mesh(data);
        set(gca, 'FontSize', fs);
        xlabel(xlabelstr, 'FontSize', labelsize);
        ylabel(ylabelstr, 'FontSize', labelsize);
        zlabel('Original data');
        moviestruct(1) = getframe(fighand);
        
        colormap('default');
        mesh(omeanmat);
        set(gca, 'FontSize', fs);
        xlabel(xlabelstr, 'FontSize', labelsize);
        ylabel(ylabelstr, 'FontSize', labelsize);
        zlabel('Overall mean');
        moviestruct(2) = getframe(fighand);
        
        
        for isvd=1:nosvd;
            SVtemp=s(isvd, isvd)*u(:, isvd)*v(:, isvd)';
            Utemp=u(:, isvd);
            Vtemp=v(:, isvd);
            
            if Vtemp<0;
                Vtemp=-1.*Vtemp;
                Utemp=-1.*Utemp;
            end;
            
            Zmax=max(max(SVtemp));
            Zmin=min(min(SVtemp));
            Vshape=Zmax*Vtemp./max(abs(Vtemp));
            Ushape=Zmax*Utemp./max(abs(Utemp));
            
            %mesh plot for the SVD components
            colormap('default');
            mesh(SVtemp);
            set(gca, 'FontSize', fs);
        xlabel(xlabelstr, 'FontSize', labelsize);
        ylabel(ylabelstr, 'FontSize', labelsize);
            zlabel(['SV' num2str(isvd)]);
            moviestruct((1+nrow+ncol)*(isvd-1)+3)=getframe(fighand);
            
            %set initial values
            x0=(1:ncol);
            y0=(1:nrow);
            yx0=zeros(1, ncol);
            xy0=zeros(1, nrow);
            yx=yx0;
            xy=xy0;
            
            for irow=1:nrow;
                yx=yx+1;
                
                colormap('gray');
                mesh(SVtemp);
                set(gca, 'FontSize', fs);
                xlabel(xlabelstr, 'FontSize', labelsize);
                ylabel([ylabelstr, ' ', num2str(irow)], 'FontSize', labelsize);
                zlabel(['SV', num2str(isvd)]);
                hold on;
                plot3(x0, yx0, Vshape, 'b-.', 'LineWidth', 4);
                plot3(xy0, y0, Ushape, 'g--', 'LineWidth', 4);
                plot3(x0, yx, SVtemp(irow, :), 'r-', 'LineWidth', 6);
                plot3(xy0(irow), y0(irow), Ushape(irow), 'r.', 'MarkerSize', 50);
                hold off;
                
                moviestruct((1+nrow+ncol)*(isvd-1)+3+irow) = getframe(fighand);
            end;
            
            for icol=1:ncol;
                xy=xy+1;
                
                colormap('gray');
                mesh(SVtemp);
                set(gca, 'FontSize', fs);
                xlabel([xlabelstr, ' ', num2str(icol)], 'FontSize', labelsize);
                ylabel(ylabelstr, 'FontSize', labelsize);
                zlabel(['SV', num2str(isvd)]);
                hold on;
                plot3(x0, yx0, Vshape, 'b-.', 'LineWidth', 4);
                plot3(xy0, y0, Ushape, 'g--', 'LineWidth', 4);
                plot3(xy, y0, SVtemp(:, icol), 'r-', 'LineWidth', 6);
                plot3(x0(icol), yx0(icol), Vshape(icol), 'r.', 'MarkerSize', 50);
                hold off;
                
                moviestruct((1+nrow+ncol)*(isvd-1)+3+nrow+icol) = getframe(fighand);
            end;
        end;

        fullorder=[(1:nrow), (nrow-1:-1:1), (1+nrow:nrow+ncol), ((nrow+ncol-1):-1:1+nrow)];
            
        initone=ones(1, 5);
        clear isvd;
        vmorder=[initone, 2*initone];
        
        for isvd=1:nosvd;
            firstposition=(1+nrow+ncol)*(isvd-1)+3;
            innerone=firstposition*initone;
            innerposition=firstposition+fullorder;
            temp=vmorder;
            vmorder=[temp innerone innerposition];
        end;
      
        moviestruct = moviestruct(vmorder);
            
        savestr='SVDmovie';
            
        movie2avi(moviestruct,savestr,'compression',moviecstr, ...
                      'keyframe',moviefps,'fps',moviefps) ;

        
        
    
    
    end; %end irow and icol
    
    
    
   
end;%end imovie