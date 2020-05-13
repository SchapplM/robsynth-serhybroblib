% Jacobian of implicit kinematic constraints of
% fourbar1turnIC
% projection from active to passive joints coordinates
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
% B21 [(no of passive joints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = fourbar1turnIC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:12
% EndTime: 2020-05-07 11:33:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (12->4), mult. (6->10), div. (7->3), fcn. (8->3), ass. (0->8)
t28 = -qJ(4) + qJ(2);
t22 = sin(qJ(3) + t28);
t21 = 0.1e1 / t22;
t27 = pkin(2) * t21 / pkin(4);
t26 = 0.1e1 / pkin(3);
t24 = sin(qJ(3));
t23 = sin(t28);
t1 = [0, (pkin(2) * t23 - pkin(3) * t22) * t26 * t21; 0, t24 * t27; 0, (t24 * pkin(3) - pkin(4) * t23) * t26 * t27;];
B21 = t1;
