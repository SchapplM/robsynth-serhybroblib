% Jacobian time derivative of explicit kinematic constraints of
% palh3m2IC
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
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_p = palh3m2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJD has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:55:19
% EndTime: 2020-05-07 04:55:19
% DurationCPUTime: 0.03s
% Computational Cost: add. (30->12), mult. (24->15), div. (0->0), fcn. (12->10), ass. (0->12)
t35 = pkin(2) * (qJD(7) + qJD(2));
t34 = pkin(7) * (qJD(7) + qJD(8));
t33 = pkin(9) * (qJD(3) + qJD(4));
t32 = pkin(3) * qJD(7);
t31 = pkin(5) * qJD(6);
t29 = -qJ(7) + pkin(15);
t25 = -qJ(8) + t29;
t30 = cos(t25) * t34;
t24 = qJ(7) + pkin(16) + qJ(2);
t23 = qJ(3) + qJ(4) + pkin(14);
t21 = sin(t25) * t34;
t1 = [0, -cos(qJ(6)) * t31, cos(t24) * t35, 0, 0, 0; 0, -sin(qJ(6)) * t31, sin(t24) * t35, 0, 0, 0; cos(t23) * t33, 0, -t30 + cos(t29) * t32, -t30, 0, 0; sin(t23) * t33, 0, t21 - sin(t29) * t32, t21, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
PhiD_p = t1;
