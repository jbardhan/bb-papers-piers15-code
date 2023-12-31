function surfdata = makeSternSrf(dielBndy, sternBndy)

surfdata = struct('dielBndy', dielBndy, 'sternBndy', sternBndy, ...
		  'numDielPanels', [length(dielBndy(1).areas)],...
		  'numSternPanels',[length(sternBndy(1).areas)]);