%{
    Copyright 2012 Shay Ben Elazar ©
    This file is part of INSP3CT.

    INSP3CT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    INSP3CT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with INSP3CT.  If not, see <http://www.gnu.org/licenses/>.
%}
% Special thanks to Florian Wagner for his implementation of mHG for Matlab this code is based on.
% This is a general implementation for the partition limited Minimum Hypergeometric test which takes a binary vector and maximum threshold and returns the corrected minimum hypergeometric p-value.
function [index,pval] = mHGthresh(v,maxthresh)
	% use this function
	% takes binary vector
	% returns cutoff index and p-value 

	[index,min_p] = mHGT(v,maxthresh);
	N = length(v);
	K = sum(v~=0);
	
	mat = zeros(K,N-K);
	
	pval = mHG_pvalue(N,K,min_p,mat);

end

function [pval] = mHG_pvalue(N,K,mHGT,mat)

	if mHGT >= 1
		pval = 1;
		return
	elseif K == 0 || K >= N
		pval = 0;
		return
	end
	
	W = N-K;
	R_ZONE = 0;
	baseHG = 1;
	mat(1,1) = 1;
	
	for n=1:(N-1)
		
		if K >= n
			min_nK = n;
			baseHG = baseHG * (K-n+1) / (N-n+1);
		else
			min_nK = K;
			baseHG = baseHG * n / (n-K);
		end
		
		tailHG = baseHG;
		currHG = baseHG;
		k = min_nK;
		w = n - min_nK;
		
		while tailHG <= mHGT+1e-16 && k > 0 && w < W
			
% 			[n,tailHG,mHGT]
			
			mat(k+1,w+1) = R_ZONE;
			currHG = currHG * (k*(N-K-n+k)) / ((n-k+1)*(K-k+1));
			tailHG = tailHG + currHG;
			
			w = w+1;
			k = k-1;
			
		end
			
		while k >= 0 && w <= W
			
			if w > 0 && k > 0
				mat(k+1,w+1) = mat(k+1,w)*(W-w+1)/(N-n+1) + mat(k,w+1)*(K-k+1)/(N-n+1);
			elseif w > 0
				mat(k+1,w+1) = mat(k+1,w)*(W-w+1)/(N-n+1);
			elseif k > 0
				mat(k+1,w+1) = mat(k,w+1)*(K-k+1)/(N-n+1);
			end
			
			w = w+1;
			k = k-1;
			
		end
		
	end
	pval = 1 - (mat(K+1,W) + mat(K,W+1));
end

function [tail] = HGT(currHG,n,N,K,k)

	min_nK = min(n,K);
	tail = currHG;
	for i=k:(min_nK-1)
		currHG = currHG*((n-1)*(K-i))/((i+1)*(N-n-K+i+1));
		tail = tail+currHG;
	end

end

function [index,mHGT] = mHGT(v,maxthresh)

	N = length(v);
	K = sum(v~=0);
	currHG = 1;
	mHGT = 1.1;
	index = 0;
	k = 0;
	
	if or (K == 0,K == N)
		index = 0;
		mHGT = 1;
		return
	end
		
	for n = 0:(maxthresh-1)
		if v(n+1) == 0
			currHG = currHG*((n+1)*(N-K-n+k))/((N-n)*(n-k+1));
			currHGT = HGT(currHG,n+1,N,K,k);
		else
			currHG = currHG*((n+1)*(K-k))/((N-n)*(k+1));
			k = k+1;
			currHGT = HGT(currHG,n+1,N,K,k);
			
			if currHGT < mHGT
				index = n;
				mHGT = currHGT;
			end
		end
	end
	
	index = index+1;

end