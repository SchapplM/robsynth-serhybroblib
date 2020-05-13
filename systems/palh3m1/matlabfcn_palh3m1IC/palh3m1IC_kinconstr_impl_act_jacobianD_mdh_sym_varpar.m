% Jacobian time derivative of explicit kinematic constraints of
% palh3m1IC
% with respect to active joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% PhiD_a [(no of constraints)x(no. of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_a = palh3m1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJD has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:27:33
% EndTime: 2020-04-20 17:27:33
% DurationCPUTime: 0.02s
% Computational Cost: add. (16->8), mult. (16->12), div. (0->0), fcn. (8->8), ass. (0->7)
t14 = pkin(2) * (qJD(7) + qJD(2));
t13 = pkin(9) * (qJD(3) + qJD(4));
t12 = pkin(1) * qJD(2);
t11 = pkin(4) * qJD(3);
t8 = pkin(16) + qJ(7) + qJ(2);
t7 = qJ(3) + qJ(4) + pkin(14);
t1 = [0, -cos(t8) * t14 - cos(qJ(2)) * t12, 0, 0; 0, -sin(t8) * t14 - sin(qJ(2)) * t12, 0, 0; 0, 0, -cos(t7) * t13 - cos(qJ(3)) * t11, 0; 0, 0, -sin(t7) * t13 - sin(qJ(3)) * t11, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
PhiD_a = t1;
