function bool = query_for_image_subset_skip(q, t, z)

bool = false;

if isempty(q)
    return;
end

if t < q.t_do_range(1) || t > q.t_do_range(2)
    bool = true; return;
elseif t >= q.t_skip_range(1) && t <= q.t_skip_range(2)
    bool = true; return;
elseif z < q.z_do_range(1) || z > q.z_do_range(2)
    bool = true; return;
elseif z >= q.z_skip_range(1) && z <= q.z_skip_range(2)
    bool = true; return;
end            