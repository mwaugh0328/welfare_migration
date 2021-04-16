function move_vec = make_movevec(movepolicy, position, params)

move_vec = [];
count = 0;

A = [];
rural_A = [eye(10), eye(10); zeros(10), zeros(10)];
urban_A = zeros(10);

for xxx = 1:params.ntypes
    
    if isequal(position,xxx)
    
        foo = squeeze(movepolicy(xxx).rural_not(15,:,:));
  
        foo = [foo(:,1),  diff(foo,1,2)];
        foo = foo(:,1:2);

        move_vec = [move_vec ; foo(:)];
        A = blkdiag(A,rural_A);
    
        foo = squeeze(movepolicy(xxx).rural_exp(79,:,:));
    
        foo = [foo(:,1),  diff(foo,1,2)];
        foo = foo(:,1:2);

        move_vec = [move_vec ; foo(:)];
        A = blkdiag(A,rural_A);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        foo = squeeze(movepolicy(xxx).urban_new(86,:,:));

        foo = [foo(:,1),  diff(foo,1,2)];
        move_vec = [move_vec ; foo(:,1)];
        A = blkdiag(A,urban_A);
    
        foo = squeeze(movepolicy(xxx).urban_old(2,:,:));
    
        foo = [foo(:,1),  diff(foo,1,2)];
        move_vec = [move_vec ; foo(:,1)];
        A = blkdiag(A,urban_A);
    end
    
end