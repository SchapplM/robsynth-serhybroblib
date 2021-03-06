% Jacobian time derivative of explicit kinematic constraints of
% palh3m1IC
% with respect to passive joint coordinates
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
% PhiD_p [(no of constraints)x(no. of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_p = palh3m1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJD has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:27:33
% EndTime: 2020-04-20 17:27:33
% DurationCPUTime: 0.03s
% Computational Cost: add. (33->14), mult. (24->15), div. (0->0), fcn. (12->10), ass. (0->12)
t37 = pkin(2) * (qJD(7) + qJD(2));
t36 = pkin(7) * (-qJD(7) - qJD(8));
t35 = pkin(9) * (qJD(3) + qJD(4));
t34 = pkin(3) * qJD(7);
t33 = pkin(5) * qJD(6);
t31 = -qJ(7) + pkin(15);
t26 = -qJ(8) + t31;
t32 = cos(t26) * t36;
t27 = pkin(16) + qJ(7) + qJ(2);
t25 = qJ(3) + qJ(4) + pkin(14);
t23 = sin(t26) * t36;
t1 = [0, cos(qJ(6)) * t33, -cos(t27) * t37, 0, 0, 0; 0, sin(qJ(6)) * t33, -sin(t27) * t37, 0, 0, 0; -cos(t25) * t35, 0, -t32 - cos(t31) * t34, -t32, 0, 0; -sin(t25) * t35, 0, t23 + sin(t31) * t34, t23, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
PhiD_p = t1;
