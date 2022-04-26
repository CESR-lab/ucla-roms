
label = '00000';
list = dir(['sample_extract.' label '*nc']);

gname = 'grid1';

nsub = length(list);
% 
% figure out files with variables for this grid
if 0
 ib = 0;
 for i = 1:nsub
   ncinf = ncinfo(list(i).name);
   fvars = ncinf.Variables;
   if length(fvars)>0
     for iv= 1:length(fvars)
       if strfind(fvars(iv).Name,gname)
         ib = ib+1;
         bry(ib).name = list(i).name;
	 break
       end
     end
   end
 end
else
 bry = list;
end

% which vars, in what files
got_time = 0;
nsub = length(bry);
iw = 0;
nvars = 0;
lvars = struct;
for i = 1:nsub

   ncinf = ncinfo(bry(i).name);
   fvars = ncinf.Variables;
   for iv= 1:length(fvars)
     if strfind(fvars(iv).Name,gname) % only if this var is part of the child grid

       if strfind(fvars(iv).Name,'time') % we only need to read one time variable
         if ~got_time                    % associated with this grid since they're 
           tname = fvars(iv).Name;       % assumed to all be the same
	   nrec = fvars(iv).Size;        % number of time records in this series of files
           tfile = bry(i).name;
           got_time = 1;
	 end
         continue                        
       end                               
       % check if we have this variable already 
       % in the list of vars
       found = 0;
       for ivar = 1:nvars
          if strcmp(lvars(ivar).name,fvars(iv).Name)
            found=1; vidx = ivar; break
          end
       end
       if ~found
         nvars = nvars+1;
         lvars(nvars).name = fvars(iv).Name;
         start = ncreadatt(bry(i).name,fvars(iv).Name,'start');
         count = ncreadatt(bry(i).name,fvars(iv).Name,'count');
         lvars(nvars).files = [i];
         lvars(nvars).start = [start];
         lvars(nvars).count = [count];
  
         dname = ncreadatt(bry(i).name,fvars(iv).Name,'dname');
         dsize = ncreadatt(bry(i).name,fvars(iv).Name,'dsize');
         lvars(nvars).dname = dname;
         lvars(nvars).dsize = dsize;

         sze = size(ncread(bry(i).name,fvars(iv).Name));
         if length(sze)>2
  	   lvars(nvars).TD = 1;
	   nz = sze(2);     % s_rho dimension size
         else
           lvars(nvars).TD = 0;
         end

       else
         start = ncreadatt(bry(i).name,fvars(iv).Name,'start');
         count = ncreadatt(bry(i).name,fvars(iv).Name,'count');
         lvars(vidx).start = [lvars(vidx).start start];
         lvars(vidx).count = [lvars(vidx).count count];
         lvars(vidx).files = [lvars(vidx).files i];
       end
     end % if name contains the gridname

   end
end

bryname = [gname '_bry.' label '.nc'] 
delete(bryname)

% create the boundary file
nccreate(bryname,'ocean_time','dimensions',{'bry_time',0},'datatype','double');
for i = 1:length(lvars)
  str = lvars(i).name;
  idx = strfind(str,'_');
  vname = [str(idx(2)+1:end) '_' str(idx(1)+1:idx(2)-1)];
  lvars(i).vname = vname;
  dname = lvars(i).dname;
  dsize = lvars(i).dsize;

  if lvars(i).TD
    nccreate(bryname,vname,'dimensions',{dname,dsize,'s_rho',nz,'bry_time',0},'datatype','single');
  else
    nccreate(bryname,vname,'dimensions',{dname,dsize,'bry_time',0},'datatype','single');
  end
end

% read the data from the partitioned files and write to the boundary file
for irec = 1:nrec

  disp(['Merging record: ',num2str(irec)])

  time = ncread(tfile,tname,[irec],[1]);
  ncwrite(bryname,'ocean_time',time,[irec]);

  for i = 1:length(lvars)
    name = lvars(i).name;
    vname = lvars(i).vname;
    if lvars(i).TD
      var = zeros(lvars(i).dsize,nz);
    else
      var = zeros(lvars(i).dsize,1);
    end

    for j = 1:length(lvars(i).files)

      idx = lvars(i).files(j);
      fname = bry(idx).name;
      count = lvars(i).count(j);
      start = lvars(i).start(j);

      if lvars(i).TD
        chunk = ncread(fname,name,[1 1 irec],[inf inf 1]);
        var(start:start+count-1,:) = chunk;
      else
        chunk = ncread(fname,name,[1 irec],[inf 1]);
        var(start:start+count-1) = chunk;
      end
    end

    if lvars(i).TD
      ncwrite(bryname,vname,var,[1 1 irec]);
    else
      ncwrite(bryname,vname,var,[1 irec]);
    end

  end
end









