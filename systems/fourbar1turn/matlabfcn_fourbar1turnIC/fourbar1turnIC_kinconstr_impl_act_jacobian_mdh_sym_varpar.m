% Jacobian of implicit kinematic constraints of
% fourbar1turnIC
% with respect to active joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% Phi_a [(no of constraints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_a = fourbar1turnIC_kinconstr_impl_act_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_act_jacobian_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_act_jacobian_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:11
% EndTime: 2020-05-07 11:33:11
% DurationCPUTime: 0.01s
% Computational Cost: add. (4->3), mult. (4->4), div. (0->0), fcn. (4->4), ass. (0->2)
t2 = qJ(2) + qJ(3);
t1 = [0, -pkin(3) * sin(t2) + sin(qJ(2)) * pkin(2); 0, pkin(3) * cos(t2) - cos(qJ(2)) * pkin(2); 0, -1;];
Phi_a = t1;
