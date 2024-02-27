function test = return_test_index(T, n_blk, blk)

% ========================================================================
% EXTRACT INDEX OF BLOCKS TO BE REMOVED FROM ESTIMATION
% >> these observations will be used for "test" sample when calculating
%    predictive densities
% ========================================================================
%
% -- INPUT --
%   T: total number of periods in sample
%   n_blocks: total number of blocks desired
%   block: iteration of which block we are omitting
%
% -- OUTPUT --
%   dropped_vec = The indices that would index data in order to produce
%
% ---
% NOTE: If block_iter = meta.n_blocks+1, we are estimating global model
% ------------------------------------------------------------------------

if blk > n_blk + 1
    error("block_iter is larger than n_blocks")
end

% NUMBER OF PERIODS IN EACH BLOCK
T_b = floor(T/n_blk); 

% OBSERVATIONS TO BE DROPPED FROM ESTIMATION
% >> USED AS "TEST SAMPLE FOR PREDICTIVE DENSITY CALCULATIONS
if blk == n_blk + 1     % global estimation
    test = 1:T;
elseif blk == n_blk     % last block
    test = (blk*T_b)-T_b+1:T;  
else                    % first (n_blocks-1) blocks
    test = (blk*T_b)-T_b+1:blk*T_b;
end
    
end
