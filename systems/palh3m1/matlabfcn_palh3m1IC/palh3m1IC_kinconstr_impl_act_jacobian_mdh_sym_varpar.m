% Jacobian of implicit kinematic constraints of
% palh3m1IC
% with respect to active joint coordinates
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
% Phi_a [(no of constraints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_a = palh3m1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_kinconstr_impl_act_jacobian_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:27:33
% EndTime: 2020-04-20 17:27:33
% DurationCPUTime: 0.02s
% Computational Cost: add. (12->6), mult. (8->8), div. (0->0), fcn. (8->8), ass. (0->3)
t6 = pkin(16) + qJ(7) + qJ(2);
t5 = qJ(3) + qJ(4) + pkin(14);
t1 = [0, -pkin(2) * sin(t6) - sin(qJ(2)) * pkin(1), 0, 0; 0, pkin(2) * cos(t6) + cos(qJ(2)) * pkin(1), 0, 0; 0, 0, -pkin(9) * sin(t5) - sin(qJ(3)) * pkin(4), 0; 0, 0, pkin(9) * cos(t5) + cos(qJ(3)) * pkin(4), 0; 0, -1, 0, 0; 0, 0, 1, 0;];
Phi_a = t1;
