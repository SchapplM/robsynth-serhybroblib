% Jacobian of implicit kinematic constraints of
% fourbar1IC
% projection from active to passive joints coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% B21 [(no of passive joints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:15
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = fourbar1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:15:48
% EndTime: 2020-04-24 20:15:48
% DurationCPUTime: 0.02s
% Computational Cost: add. (12->4), mult. (6->10), div. (7->3), fcn. (8->3), ass. (0->8)
t28 = -qJ(3) + qJ(1);
t22 = sin(qJ(2) + t28);
t21 = 0.1e1 / t22;
t27 = pkin(2) * t21 / pkin(4);
t26 = 0.1e1 / pkin(3);
t24 = sin(qJ(2));
t23 = sin(t28);
t1 = [(pkin(2) * t23 - pkin(3) * t22) * t26 * t21; t24 * t27; (t24 * pkin(3) - pkin(4) * t23) * t26 * t27;];
B21 = t1;
