% Jacobian of implicit kinematic constraints of
% palh3m2IC
% with respect to passive joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% Phi_p [(no of constraints)x(no of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_p = palh3m2IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:55:19
% EndTime: 2020-05-07 04:55:19
% DurationCPUTime: 0.02s
% Computational Cost: add. (23->9), mult. (12->10), div. (0->0), fcn. (12->10), ass. (0->7)
t20 = -qJ(7) + pkin(15);
t19 = -qJ(8) + t20;
t18 = qJ(7) + pkin(16) + qJ(2);
t17 = qJ(3) + qJ(4) + pkin(14);
t16 = cos(t19) * pkin(7);
t15 = sin(t19) * pkin(7);
t1 = [0, -sin(qJ(6)) * pkin(5), pkin(2) * sin(t18), 0, 0, 0; 0, cos(qJ(6)) * pkin(5), -pkin(2) * cos(t18), 0, 0, 0; pkin(9) * sin(t17), 0, t15 - pkin(3) * sin(t20), t15, 0, 0; -pkin(9) * cos(t17), 0, t16 - pkin(3) * cos(t20), t16, 0, 0; 0, 1, -1, 0, -1, 0; 1, 0, -1, -1, 0, 1;];
Phi_p = t1;
