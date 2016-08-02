select
	substring(hh.hmm_acc,1,12), substring(hh.hmm_com_name,1,45), er.role_id, substring(er.mainrole,1,35), substring(er.sub1role,1,45), hh.ec_num
from 
	hmm..hmm3 hh, hmm..hmm_role_link hrl, egad..roles er
where 
	hh.hmm_acc in 
       (select  h.hmm_acc
          from prop_step ps, step_ev_link se, hmm..hmm3 h, prop_marker m
          where ps.prop_def_id = 64528
            and ps.prop_step_id *= se.prop_step_id
            and se.query *= h.hmm_acc
            and se.query *= m.query
            and m.prop_def_id = 64528
            and h.is_current = 1)
	and hh.hmm_acc = hrl.hmm_acc
	and hrl.role_id = er.role_id
order by 
	er.mainrole

