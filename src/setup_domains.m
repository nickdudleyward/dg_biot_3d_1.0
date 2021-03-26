function param = setup_domains(param)

% param.elt_to_domain holds domain codes, indexed by elements
% domain_to_model_type (extracted from param.phys.model_type) holds model 
% type codes, indexed by domain codes
% domain codes are only used as indices into this array and can be assigned
% arbitrarily
% It is essential that domain_to_model_type(elt_to_domain) evaluates to a
% row vector (this error now trapped but should be further investigated)

% model type codes are as defined in setup_magic
% TODO should this code be in build_params?

Globals3D;

% Find model types for each element

domain_to_model_type = [ param.phys.model_type ];
param.elt_to_model_type = domain_to_model_type(param.elt_to_domain);

% Everything goes horribly wrong without actual errors if this is a column 
% vector instead of a row vector, so check at this point (and some time
% work out what actually happens!)

if size(param.elt_to_model_type, 1) ~= 1
    error('domain_to_model_type(param.elt_to_domain) does not evaluate to row vector')
end

% Extract lists of element numbers associated with different model types.
% First two lines work because all elastic model type numbers are 
% < param.POROELASTIC and all poroelastic model type numbers are 
% > param.POROELASTIC. 

param.elastic_elt     = find(param.elt_to_model_type < param.POROELASTIC);
param.poroelastic_elt = find(param.elt_to_model_type > param.POROELASTIC);
param.inviscid_elt    = find(param.elt_to_model_type == param.INVISCID);
param.lf_elt          = find(param.elt_to_model_type == param.LOW_FREQUENCY);
param.hf_elt          = find(param.elt_to_model_type == param.HIGH_FREQUENCY);

% Now look at subtypes. Find the list of index numbers within
% param.poroelastic_elt corresponding to the three subtypes inviscid, low
% frequency, high frequency

poro_model_type = param.elt_to_model_type(param.poroelastic_elt);
param.in_within_poro = find(poro_model_type == param.INVISCID);
param.lf_within_poro = find(poro_model_type == param.LOW_FREQUENCY);
param.hf_within_poro = find(poro_model_type == param.HIGH_FREQUENCY);

% Now inflate elt_to_model_type to a list indexed by element node numbers.
% First duplicate the type list Np (number of element nodes per element) 
% times, as rows of a matrix (so param.elt_node_to_model_type has the same 
% indexing structure as e.g. x, y, z)

param.elt_node_to_model_type = repmat(param.elt_to_model_type,[Np, 1]);

% pe_mask and ep_mask are indexed by face node numbers. 
% pe-mask holds true for each face node with elastic element and 
%  corresponding exterior poroelastic element
% ep-mask holds true for each face node with poroelastic element and 
%  corresponding exterior elastic element

%save('vmap.mat', 'vmapM', 'vmapP');

pe_mask = (param.elt_node_to_model_type(vmapM) < param.POROELASTIC) & ...
          (param.elt_node_to_model_type(vmapP) > param.POROELASTIC);

ep_mask = (param.elt_node_to_model_type(vmapM) > param.POROELASTIC) & ...
          (param.elt_node_to_model_type(vmapP) < param.POROELASTIC);

% pe_face_node is a list of face node numbers whose element is elastic 
% and whose corresponding exterior element is poroelastic. 
% ep_face_node is a list of face node numbers whose element is poroelastic 
% and whose corresponding exterior element is elastic. 
% Note that boundary nodes have exterior and interior nodes equal so have 
% the same model type and are hence not present in this list. This is used 
% to insert flux data in the right place. Coordinates of interface nodes 
% are e.g. x(vmapM(param.ep_face_node))

param.ep_face_node = find(ep_mask);
param.pe_face_node = find(pe_mask);

% Similar thing for elastic-elastic and poroelastic-poroelastic, with
% one major difference: a pair of (elastic, poroelastic) nodes appears
% once in pe_face_node (elastic as interior, poroelastic as exterior) 
% and once in ep_face_node (poroelastic as interior, elastic as exterior)
% in ee_face_node and pp_face_node, each pair of nodes appears twice,
% with both possible interior-exterior roles

ee_mask = (param.elt_node_to_model_type(vmapM) < param.POROELASTIC) & ...
          (param.elt_node_to_model_type(vmapP) < param.POROELASTIC);
pp_mask = (param.elt_node_to_model_type(vmapM) >  param.POROELASTIC) & ...
          (param.elt_node_to_model_type(vmapP) >  param.POROELASTIC);

param.ee_face_node = find(ee_mask);
param.pp_face_node = find(pp_mask);

% Note that pe_face_node and ep_face_node do not list corresponding points
% in the same order: pe_face_node[n] and ep_face_node[n] do not usually
% have the same location in space (although their faces do).

% The following alternative construction of ep_face_node gives rise to
% the same face node numbers, but in a different order, such that 
% pe_face_node[n] and ep_face_node[n] have the same location in space
% (as do their faces).

% if param.align_ep_pe_nodes is true, the alternative version replaces the
% value of ep_face_node calcualted above

% This order is required by get_face_flux_ep_pe, which constructs face 
% fluxes both into and out of an elastic element, only data derived from 
% that element. get_face_flux_ep, in contrast, can use either node ordering

if param.align_ep_pe_nodes

    % Construct mapping from each face node to its corresponding face node
    % on adjacent face. Boundary node numbers map to themselves.

    % node_code has one row per face node and two columns.
    % In each row, the first column of contains a face node number and
    % the second column contains an integer encoding the two element node
    % numbers associated with that face node number (or, for boundary nodes, 
    % two copies of the single face node number). The code formula
    % describes an injective map from unordered pairs of element node numbers.

    num_face_nodes = numel(vmapM);
    min_en = min([vmapM, vmapP],[],2);
    max_en = max([vmapM, vmapP],[],2);
    node_code = [ (1:num_face_nodes)', max_en.*(max_en-1)/2+min_en-1 ];

    % Sort the rows of node_node by code number

    node_code = sortrows(node_code, 2);

    % Corresponding face node numbers (column 1) now appear in consecutive 
    % rows with the same value in column 2 (because they have the same vmapP 
    % and vmapM values, but in different orders; the node code in column 2 is 
    % insensitive to order). This allows us to build a mapping from face node 
    % number to corresponding exterior face node number.
    % Boundary nodes need to be handled separately: they appear once only
    % and the code number in column 2 is unique.

    param.face_node_map = zeros(1,num_face_nodes);
    j=1;
    while j<=num_face_nodes
        if j<num_face_nodes && node_code(j,2)==node_code(j+1,2)
            param.face_node_map(node_code(j,1))   = node_code(j+1,1);
            param.face_node_map(node_code(j+1,1)) = node_code(j,1);
            j = j+1;
        else
            param.face_node_map(node_code(j,1)) = node_code(j,1);
        end
        j = j+1;
    end

    % Quick sanity check: each face node number should appear exactly once

    if any(sort(param.face_node_map) ~= 1:num_face_nodes)
        error('internal error: face_node_map is not a bijection')
    end

    % needs to be a column vector

    param.ep_face_node = param.face_node_map(param.pe_face_node)';

end % if param.align_ep_pe_nodes
