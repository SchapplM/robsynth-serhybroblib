% Jacobian of implicit kinematic constraints of
% picker2Dm2IC
% with respect to active joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% Phi_a [(no of constraints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 09:21
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_a = picker2Dm2IC_kinconstr_impl_act_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2IC_kinconstr_impl_act_jacobian_mdh_sym_varpar: qJ has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2IC_kinconstr_impl_act_jacobian_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 09:20:35
% EndTime: 2020-05-11 09:20:35
% DurationCPUTime: 0.02s
% Computational Cost: add. (14->9), mult. (12->10), div. (0->0), fcn. (12->10), ass. (0->6)
t15 = sin(qJ(1)) * pkin(1);
t14 = cos(qJ(1)) * pkin(1);
t11 = qJ(1) + qJ(2);
t10 = qJ(1) + qJ(8);
t9 = qJ(4) + t11;
t1 = [-pkin(4) * sin(t9) - pkin(3) * sin(t11) - t15, cos(qJ(7)) * pkin(3); pkin(4) * cos(t9) + pkin(3) * cos(t11) + t14, sin(qJ(7)) * pkin(3); -pkin(5) * sin(t10) + t15, 0; pkin(5) * cos(t10) - t14, 0; 0, 0; 0, 0; 1, -1; -1, 0; 0, 0; 0, 0;];
Phi_a = t1;
