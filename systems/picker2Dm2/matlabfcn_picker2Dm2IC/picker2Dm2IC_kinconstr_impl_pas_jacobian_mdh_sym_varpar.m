% Jacobian of implicit kinematic constraints of
% picker2Dm2IC
% with respect to passive joint coordinates
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
% Phi_p [(no of constraints)x(no of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 09:21
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_p = picker2Dm2IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: qJ has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 09:20:35
% EndTime: 2020-05-11 09:20:35
% DurationCPUTime: 0.08s
% Computational Cost: add. (27->14), mult. (18->14), div. (0->0), fcn. (18->14), ass. (0->10)
t39 = qJ(1) + qJ(2);
t35 = qJ(4) + t39;
t41 = pkin(4) * sin(t35);
t37 = qJ(3) + qJ(9);
t40 = cos(t37) * pkin(2);
t38 = qJ(1) + qJ(8);
t36 = pkin(8) + qJ(5);
t32 = sin(t37) * pkin(2);
t31 = pkin(4) * cos(t35);
t1 = [-t41 - pkin(3) * sin(t39), 0, -t41, 0, 0, 0, 0, 0, 0, 0; t31 + pkin(3) * cos(t39), 0, t31, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, sin(t36) * pkin(1), 0, -pkin(5) * sin(t38), 0, 0, 0, 0; 0, 0, 0, -cos(t36) * pkin(1), 0, pkin(5) * cos(t38), 0, 0, 0, 0; 0, t32 - sin(qJ(3)) * pkin(6), 0, 0, sin(qJ(6)) * pkin(6), 0, t32, 0, 0, 0; 0, -t40 + cos(qJ(3)) * pkin(6), 0, 0, -cos(qJ(6)) * pkin(6), 0, -t40, 0, 0, 0; 1, 0, 1, 0, 0, 0, 0, 1, 0, 0; 0, 0, 0, 1, 0, -1, 0, 0, -1, 0; 0, -1, 0, 0, 1, 0, -1, 0, 0, 1; 1, 0, 0, 0, 0, -1, -1, 0, 0, 0;];
Phi_p = t1;
