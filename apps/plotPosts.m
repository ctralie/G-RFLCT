function [] = plotPosts(posts)
	posts = posts';
	pos = [-1, 0, 1];
	N = size(posts, 2);
	for i = 1:3
		for j = 1:3
			%fprintf(1, 'i = %i, j = %i', i, j)
			k = (i-1)*3+j;
			subplot(3, 3, k);
			plot(posts(k, :));
			axis([1, N, 0, 1]);
			title(sprintf('Source Position (%i, %i, -2)', pos(i), pos(j)));
		end
	end
	lastPost = posts(:, N);
	[C, Index] = max(lastPost);
	fprintf(1, 'Strongest is %i\n', Index-1)
end
