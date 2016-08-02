select 
	s.step_num, e.query, g.go_term
from 
	prop_step s, step_ev_link e, hmm..hmm_go_link g 
where 
	prop_def_id = 64528 and 
	s.prop_step_id = e.prop_step_id and 
	e.query = g.hmm_acc 
		order by g.hmm_acc, g.go_term

